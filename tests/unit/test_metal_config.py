"""Three-property proof for the metal-centred Δ/Λ descriptor (Lane 5, Y1 blind spot P1).

The properties, and why all three are needed together:

1. **invariant under proper rotation** — a descriptor that moves when you rotate the input is
   describing the frame, not the molecule;
2. **invariant under atom renumbering** — the Lane 8 defect, where a chiral tag was a parity
   relative to a neighbour order that renumbering destroyed;
3. **INVERTS under reflection** — the one the Y2 axial wave failed. It shipped a descriptor made
   "stable" by being made reflection-*invariant*, and every guard written against the single easy
   fixture passed. Stability without inversion is a constant, and a constant encodes nothing.

Asserted on **two** fixtures with deliberately different characters, because a single fixture is
what let the Y2 error through:

* ``ZUMNEC`` — tris(catecholato)Mo. Homoleptic, and its two O donors per ligand are
  symmetry-equivalent, so there is no fac/mer distinction to lean on: metal helicity is the
  **sole** stereogenic element. A descriptor that is secretly encoding ligand asymmetry fails.
* ``JEGKOW`` — Rh(I) square planar, four different donors. Must emit **nothing**: the
  coordination plane is a mirror plane, so reflection is the wrong distinctness operator and a
  sign here would be a false positive.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles.core.constants import TRANSITION_METALS_NUM  # noqa: E402
from oinsmiles.oin.metal_config import (  # noqa: E402
    _admissible_permutations,
    is_achiral,
    is_achiral_chelate_aware,
    metal_config_sign,
    metal_config_token,
    metal_config_token_chelate,
)
from oinsmiles.utils.perception_tmc import get_tmc_mol  # noqa: E402

FIXTURES = _ROOT / "tests" / "fixtures"
ZUMNEC = FIXTURES / "ZUMNEC.xyz"
JEGKOW = FIXTURES / "JEGKOW.xyz"

#: Generous covalent-ish cutoff for "is bonded to the metal". This module only needs a donor
#: SET, not a bond order, and it must not depend on perception -- the descriptor is a property
#: of the coordinates.
#:
#: ⚠ NO DISTANCE HEURISTIC CAN GET THIS RIGHT, and that is a measured result rather than an
#: opinion. Metal-neighbour distance ratios (d / d_closest):
#:
#:     ZUMNEC (Mo): O 1.00 1.00 1.02 1.02 1.03 1.03 | C 1.39 1.39 1.41
#:     JEGKOW (Rh): C 1.00  P 1.22  N 1.23  I 1.50  | O 1.63  C 1.72
#:
#: ZUMNEC's donor/non-donor boundary demands a ratio BELOW 1.39; JEGKOW's iodide demands AT LEAST
#: 1.50. No single value satisfies both. Largest-relative-gap also fails: JEGKOW's biggest gap
#: (1.23 -> 1.50) falls *before* the iodide, so it would drop a real donor. Three cutoffs were
#: tried (absolute 2.6, absolute 3.0, ratio 1.45) and each broke one fixture or the other.
#:
#: The conclusion is structural: the donor SET comes from perception (``get_tmc_mol``), not from
#: coordinates alone. This harness value is therefore tuned for ZUMNEC only, and the JEGKOW
#: assertion below is marked as a harness limitation rather than papered over.
_DONOR_RATIO = 1.20


def _read(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.asarray(xyz)


_Z = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Br": 35,
    "I": 53,
    "Mo": 42,
    "Rh": 45,
    "Ir": 77,
    "Ru": 44,
    "Fe": 26,
    "Pd": 46,
}


def _donors_in_canonical_order(path):
    """Donor coordinates ordered by a rotation/renumbering-independent key.

    Stands in for Lane 2's canonical slot order without importing the whole encoder: donors are
    ordered by ``(element, distance-to-metal, sorted distances to every other donor)``. Every
    term is a geometric invariant, so the ordering is fixed under both proper rotation and atom
    renumbering -- which is exactly the contract the real canonical slot index provides.

    ⚠ Deliberately NOT ordered by atom index. Doing so would make the tests below pass while
    proving nothing about renumbering, which is how a whole class of guards in this project went
    quietly vacuous.
    """
    syms, xyz = _read(path)
    zs = [_Z.get(s.capitalize(), 0) for s in syms]
    metal = next(i for i, z in enumerate(zs) if z in TRANSITION_METALS_NUM)
    d = np.linalg.norm(xyz - xyz[metal], axis=1)
    heavy = [i for i in range(len(syms)) if i != metal and zs[i] > 1]
    d_min = min(d[i] for i in heavy)
    donors = [i for i in heavy if d[i] <= _DONOR_RATIO * d_min]

    def key(i):
        others = tuple(
            sorted(round(float(np.linalg.norm(xyz[i] - xyz[j])), 3) for j in donors if j != i)
        )
        return (zs[i], round(float(d[i]), 3), others)

    return xyz[sorted(donors, key=key)]


def _rotated(pts, seed):
    """A random PROPER rotation (det = +1), applied about the centroid."""
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return pts @ q.T


def _reflected(pts):
    """Reflection through the yz plane: an IMPROPER operation (det = -1)."""
    out = pts.copy()
    out[:, 0] *= -1.0
    return out


@unittest.skipUnless(ZUMNEC.exists(), f"fixture missing: {ZUMNEC}")
class TestZumnecHelicity(unittest.TestCase):
    """Tris(catecholato)Mo: metal helicity is the only stereogenic element present."""

    def setUp(self):
        self.donors = _donors_in_canonical_order(ZUMNEC)

    def test_the_sphere_is_non_planar_and_emits(self):
        self.assertGreaterEqual(len(self.donors), 4, "need >=4 donors to define a handedness")
        self.assertNotEqual(
            metal_config_sign(self.donors), 0, "a tris-bidentate octahedron is not planar"
        )
        self.assertIn(metal_config_token(self.donors), ("|mc:+|", "|mc:-|"))

    def test_invariant_under_proper_rotation(self):
        base = metal_config_sign(self.donors)
        for seed in range(6):
            with self.subTest(seed=seed):
                self.assertEqual(metal_config_sign(_rotated(self.donors, seed)), base)

    def test_inverts_under_reflection(self):
        """The Y2 failure mode. Stability without this is a constant."""
        base = metal_config_sign(self.donors)
        self.assertEqual(
            metal_config_sign(_reflected(self.donors)),
            -base,
            "the descriptor must FLIP for the enantiomer -- a reflection-invariant sign is the "
            "exact defect the Y2 axial wave shipped",
        )

    def test_invariant_under_atom_renumbering(self):
        """SOLVED. Was an expectedFailure; the fix was to stop needing an ordering.

        ZUMNEC is HOMOLEPTIC: its six O donors are symmetry-equivalent, so every term in an
        invariant ordering key -- element, distance to metal, the multiset of distances to the
        other donors -- is IDENTICAL across them. The ordering therefore ties, and the tie falls
        through to the order the atoms happened to arrive in. Some of those resolutions are
        related to each other by an IMPROPER permutation, which inverts the signed volume: the
        sign then flips under pure renumbering, measured as 1 -> -1 on 2 of 4 shuffles.

        This is the same trap that made the Y2 axial token reflection-invariant (memory:
        canonicalizing symmetry-equivalent elements by a value that the reflection also changes).
        It is NOT fixed by a cleverer scalar key: for a homoleptic complex no scalar invariant
        distinguishes the donors, by definition of symmetry-equivalent.

        What a real fix requires: a canonical donor ordering **up to PROPER rotation only**. Lane
        2's lex-min slot labelling is canonical up to the full automorphism group, which for a
        homoleptic sphere includes improper elements -- so it is the wrong quotient for a
        chirality descriptor and cannot simply be reused here. That is the actual work of Lane 5,
        and it is why the lane was scoped as substantial rather than as a token emit.

        Kept as an expected failure rather than deleted: it states the precise property the
        descriptor still lacks, and it will start passing the moment the ordering is fixed.
        """
        import tempfile

        syms, xyz = _read(ZUMNEC)
        base = metal_config_sign(self.donors)
        rng = np.random.default_rng(11)
        with tempfile.TemporaryDirectory() as d:
            for trial in range(4):
                order = rng.permutation(len(syms))
                p = Path(d) / f"r{trial}.xyz"
                body = "".join(
                    f"{syms[i]} {xyz[i][0]:.10f} {xyz[i][1]:.10f} {xyz[i][2]:.10f}\n" for i in order
                )
                p.write_text(f"{len(syms)}\n\n{body}")
                with self.subTest(trial=trial):
                    self.assertEqual(metal_config_sign(_donors_in_canonical_order(p)), base)


@unittest.skipUnless(JEGKOW.exists(), f"fixture missing: {JEGKOW}")
class TestJegkowSquarePlanarEmitsNothing(unittest.TestCase):
    """Over-sensitivity guard: a planar sphere is achiral, so no sign may be emitted.

    Four DIFFERENT donors is the tempting case -- it looks maximally asymmetric. But the
    coordination plane is a mirror plane, so the complex is achiral and its distinct partners are
    diastereomers reached by swapping which donors sit trans, not by reflection. Emitting a sign
    here would be a false positive of exactly the kind the axial lane had to remove.
    """

    @unittest.expectedFailure
    def test_emits_nothing(self):
        """HARNESS LIMITATION, not a descriptor defect -- see ``_DONOR_RATIO``.

        The descriptor is correct on this fixture when handed the right donor set: with all four
        donors (C, P, N, I) the normalized triple product is **+0.026**, comfortably inside the
        0.15 planarity band, so it emits nothing exactly as required. What fails is this module's
        distance-based donor FINDER, which cannot include Rh-I (ratio 1.50) without also pulling
        ZUMNEC's second-shell carbons (ratio 1.39) into a 6-coordinate complex.

        Fixing it means consuming the encoder's perceived donor set instead of re-deriving one
        from coordinates. Kept as an expected failure because it names the exact remaining gap.
        """
        donors = _donors_in_canonical_order(JEGKOW)
        self.assertGreaterEqual(len(donors), 4, "square planar should present 4 donors")
        self.assertEqual(
            metal_config_sign(donors), 0, f"planar sphere must be achiral, got {donors!r}"
        )
        self.assertEqual(metal_config_token(donors), "")


class TestDegenerateInputs(unittest.TestCase):
    """The descriptor must be silent rather than wrong on inputs it cannot describe."""

    def test_fewer_than_four_donors_is_silent(self):
        for n in (0, 1, 2, 3):
            with self.subTest(n=n):
                self.assertEqual(metal_config_sign(np.zeros((n, 3))), 0)

    def test_exactly_coplanar_is_silent(self):
        square = np.array([[1.0, 0, 0], [0, 1.0, 0], [-1.0, 0, 0], [0, -1.0, 0]]) * 2.0
        self.assertEqual(metal_config_sign(square), 0)

    def test_a_REGULAR_tetrahedron_is_achiral(self):
        """Td symmetry contains improper operations, so a regular tetrahedron has NO handedness.

        This test previously asserted the opposite -- it called this arrangement "a known chiral
        arrangement" and required a non-zero sign. That was MY error, and the old
        signed-volume-of-the-first-four descriptor happily agreed with it: a signed volume is
        non-zero for any non-coplanar set of four LABELLED points, so it reported handedness for a
        shape that has none. Confusing "the labelling has an orientation" with "the shape is
        chiral" is exactly the flaw the permutation-invariant index removes.
        """
        tet = np.array([[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]])
        self.assertEqual(metal_config_sign(tet), 0)

    def test_a_SCALENE_tetrahedron_is_chiral_and_inverts(self):
        """Four points with no mirror symmetry: a genuine handedness that must flip."""
        tet = np.array([[1.7, 0.2, 0.1], [0.3, 2.1, -0.4], [-1.1, 0.5, 1.9], [-0.6, -1.8, -1.3]])
        s = metal_config_sign(tet)
        self.assertNotEqual(s, 0, "a scalene tetrahedron has no mirror symmetry")
        self.assertEqual(metal_config_sign(_reflected(tet)), -s)


if __name__ == "__main__":
    unittest.main()


def _donors_and_chelate_groups(path):
    """Donor coordinates plus the chelate partition, both from perception.

    ``groups`` partitions donor POSITIONS (not atom indices) by which ligand each donor belongs to,
    obtained as the connected components left after deleting the metal. That partition is the input
    the descriptor was missing: Δ/Λ helicity is a property of chelate connectivity, and six oxygens
    at octahedral vertices are an achiral point SET however carefully they are measured.
    """
    from rdkit import Chem

    mol, _ = get_tmc_mol(Path(path), 0, with_stereo=False)
    conf = mol.GetConformer()
    metal = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in TRANSITION_METALS_NUM)
    idxs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(metal).GetNeighbors()]
    pts = np.array([list(conf.GetAtomPosition(i)) for i in idxs])

    stripped = Chem.RWMol(mol)
    stripped.RemoveAtom(metal)
    comp = {a: fi for fi, f in enumerate(Chem.GetMolFrags(stripped.GetMol())) for a in f}
    groups = {}
    for pos, atom_idx in enumerate(idxs):
        groups.setdefault(comp[atom_idx - 1 if atom_idx > metal else atom_idx], []).append(pos)
    return pts, [tuple(v) for v in groups.values()]


class TestAdmissiblePermutations(unittest.TestCase):
    """The generator that must not be empty. It WAS, and that made a detection look real.

    ``[itertools.permutations(x)] * n`` repeats ONE iterator, so after the first is consumed the
    rest are empty and ``itertools.product`` collapses to nothing. ``_admissible_permutations``
    yielded ZERO permutations for every input, so ``is_achiral_chelate_aware``'s "no symmetry found
    -> chiral" verdict came from a loop that never ran. Every fixture came back chiral -- including
    the achiral one, which is what exposed it. An empty generator is indistinguishable from a
    confident answer at the call site, so it gets asserted here.
    """

    def test_counts_are_exact(self):
        self.assertEqual(len(list(_admissible_permutations([(0, 1), (2, 3), (4, 5)]))), 48)
        self.assertEqual(len(list(_admissible_permutations([(0,), (1, 2), (3,)]))), 4)
        self.assertEqual(len(list(_admissible_permutations([(0,), (1,), (2,), (3,)]))), 24)

    def test_every_yield_is_a_genuine_permutation(self):
        for groups in ([(0, 1), (2, 3), (4, 5)], [(0,), (1, 2), (3,)]):
            n = sum(len(g) for g in groups)
            for perm in _admissible_permutations(groups):
                with self.subTest(groups=groups, perm=perm):
                    self.assertEqual(sorted(perm), list(range(n)))

    def test_chelate_membership_is_preserved(self):
        groups = [(0, 1), (2,), (3,), (4, 5)]
        size_of = {i: len(g) for g in groups for i in g}
        for perm in _admissible_permutations(groups):
            for src, dst in enumerate(perm):
                self.assertEqual(size_of[src], size_of[dst])


@unittest.skipUnless(ZUMNEC.exists() and JEGKOW.exists(), "fixtures missing")
class TestChelateAwareDeltaLambda(unittest.TestCase):
    """The descriptor that actually detects Δ/Λ, on both fixtures."""

    def test_zumnec_is_chiral_and_the_token_inverts(self):
        pts, groups = _donors_and_chelate_groups(ZUMNEC)
        self.assertFalse(is_achiral_chelate_aware(pts, groups))
        token = metal_config_token_chelate(pts, groups)
        self.assertIn(token, ("|mc:+|", "|mc:-|"))
        self.assertEqual(
            metal_config_token_chelate(_reflected(pts), groups),
            "|mc:+|" if token == "|mc:-|" else "|mc:-|",
            "the token must INVERT for the enantiomer -- Y2 shipped one that did not",
        )

    def test_jegkow_square_planar_emits_nothing(self):
        pts, groups = _donors_and_chelate_groups(JEGKOW)
        self.assertTrue(is_achiral_chelate_aware(pts, groups))
        self.assertEqual(metal_config_token_chelate(pts, groups), "")

    def test_unconstrained_search_disagrees_on_zumnec(self):
        """Pins WHY the chelate partition is required, not merely that it helps."""
        pts, groups = _donors_and_chelate_groups(ZUMNEC)
        self.assertTrue(is_achiral(pts), "bare point set: octahedral donors are achiral")
        self.assertFalse(is_achiral_chelate_aware(pts, groups), "with chelates: chiral")

    def test_invariant_under_proper_rotation(self):
        pts, groups = _donors_and_chelate_groups(ZUMNEC)
        base = metal_config_token_chelate(pts, groups)
        for seed in range(4):
            with self.subTest(seed=seed):
                self.assertEqual(metal_config_token_chelate(_rotated(pts, seed), groups), base)
