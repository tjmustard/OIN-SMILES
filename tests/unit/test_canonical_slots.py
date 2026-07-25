"""Guards for ``oin.canonical_slots`` -- the proper-rotation group and lex-min relabeling.

These tests exist because the module is load-bearing twice over: it decides the comparison
KEY's fac/mer awareness (since v0.4.4) and, from v0.4.5, the canonical ``{n}`` slot numbers
the encoder emits. A bug here either reintroduces fac/mer blindness (silently merging two
real isomers) or destroys stereochemistry by folding over a reflection.
"""

import itertools
import unittest

import numpy as np

from oinsmiles.oin.canonical_slots import (
    GEOMETRY_VERTICES,
    canonical_slot_permutation,
    derive_rotation_group,
    geometry_rotation_group,
    lexmin_vertex_signature,
)
from oinsmiles.utils.oin_aligner import TEMPLATE_SPECS, OINDiscreteAligner

# Proper-rotation group order per geometry. These are the point group's ROTATIONAL
# subgroup orders (chiral/proper only -- no reflections, no inversion):
#   LIN D_inf_h -> the 2 vertices swap under a C2 perpendicular to the axis
#   TPL D3h     -> C3v rotations acting on 3 coplanar vertices: all 3! = 6 (planar set,
#                  see derive_rotation_group's rank<3 branch)
#   SPL D4h     -> 8 (the planar set admits the full dihedral group of the square)
#   TET Td      -> T = 12
#   TPY C3v     -> C3 = 3
#   TBP D3h     -> D3 = 6
#   SPY C4v     -> C4 = 4
#   OCT Oh      -> O = 24
#   PBP D5h     -> D5 = 10
#   SQA D4d     -> C4 = 4 (the square antiprism's proper rotations that fix the vertex
#                  partition as tabulated)
#   TCT D3h     -> C3 x C2 = ... as tabulated: 2
EXPECTED_GROUP_ORDER = {
    "LIN": 2,
    "TPL": 6,
    "SPL": 8,
    "TET": 12,
    "TPY": 3,
    "TBP": 6,
    "SPY": 4,
    "OCT": 24,
    "PBP": 10,
    "SQA": 4,
    "TCT": 2,
}


class TestVertexTableIsSingleSourceOfTruth(unittest.TestCase):
    """TD-005: ``GEOMETRY_VERTICES`` and ``TEMPLATE_SPECS['pos']`` must not drift apart."""

    def test_same_geometry_keys(self):
        self.assertEqual(set(GEOMETRY_VERTICES), set(TEMPLATE_SPECS))

    def test_same_vertex_directions(self):
        for geo, verts in GEOMETRY_VERTICES.items():
            spec = TEMPLATE_SPECS[geo]
            self.assertEqual(len(verts), len(spec), f"{geo}: vertex count differs")
            for slot in sorted(spec):
                a = np.asarray(verts[slot], dtype=float)
                b = np.asarray(spec[slot]["pos"], dtype=float)
                a = a / np.linalg.norm(a)
                b = b / np.linalg.norm(b)
                np.testing.assert_allclose(
                    a, b, atol=1e-4, err_msg=f"{geo} slot {slot}: direction differs"
                )


class TestRotationGroupIsAGroup(unittest.TestCase):
    def test_expected_orders(self):
        for geo, want in EXPECTED_GROUP_ORDER.items():
            self.assertEqual(len(geometry_rotation_group(geo)), want, f"{geo}")

    def test_elements_are_permutations(self):
        for geo, verts in GEOMETRY_VERTICES.items():
            n = len(verts)
            for perm in geometry_rotation_group(geo):
                self.assertEqual(sorted(perm), list(range(n)), f"{geo}: {perm} not a permutation")

    def test_contains_identity(self):
        for geo, verts in GEOMETRY_VERTICES.items():
            self.assertIn(tuple(range(len(verts))), geometry_rotation_group(geo), geo)

    def test_closed_under_composition_and_inverses(self):
        for geo in GEOMETRY_VERTICES:
            group = set(geometry_rotation_group(geo))
            for p in group:
                inv = [0] * len(p)
                for i, img in enumerate(p):
                    inv[img] = i
                self.assertIn(tuple(inv), group, f"{geo}: {p} has no inverse in the group")
            for p, q in itertools.product(group, repeat=2):
                comp = tuple(q[p[i]] for i in range(len(p)))
                self.assertIn(comp, group, f"{geo}: {p} o {q} left the group")

    def test_spanning_geometries_admit_only_proper_rotations(self):
        """The whole point: a reflection maps a structure to its MIRROR.

        Folding over an improper operation would collapse enantiomers -- destroying exactly
        the stereochemistry v0.4.5 exists to capture -- and would flip the eta winding sense.
        For every rank-3 geometry, reconstruct each permutation's orthogonal map and require
        ``det == +1``.
        """
        for geo, verts in GEOMETRY_VERTICES.items():
            V = np.asarray(verts, dtype=float)
            V = V / np.linalg.norm(V, axis=1, keepdims=True)
            if np.linalg.matrix_rank(V, tol=1e-6) < 3:
                continue  # planar/linear: the free out-of-plane axis is handled separately
            basis = []
            for i in range(len(V)):
                if np.linalg.matrix_rank(V[basis + [i]], tol=1e-6) == len(basis) + 1:
                    basis.append(i)
                if len(basis) == 3:
                    break
            binv = np.linalg.inv(V[basis].T)
            for perm in geometry_rotation_group(geo):
                rot = V[[perm[b] for b in basis]].T @ binv
                self.assertGreater(np.linalg.det(rot), 0, f"{geo}: {perm} is an IMPROPER operation")

    def test_octahedron_is_24_not_48(self):
        """A guard with teeth: 48 would mean reflections leaked into the group."""
        self.assertEqual(len(geometry_rotation_group("OCT")), 24)

    def test_unknown_geometry_is_none(self):
        self.assertIsNone(geometry_rotation_group("NOPE"))


class TestBruteForceSymmetriesIsIncomplete(unittest.TestCase):
    """Pins the measured gap that justifies unifying on ``derive_rotation_group``.

    ``oin_aligner._brute_force_symmetries`` builds the same group by brute-forcing Euler
    angles from the fixed grid ``[0, 90, 120, 180, 240, 270]``. That grid cannot express
    PBP's 72-degree five-fold equatorial rotation, so it finds only 2 of the 10 proper
    rotations -- the encoder currently cannot canonicalize a pentagonal-bipyramidal
    equatorial labeling at all.
    """

    @staticmethod
    def _brute(geo):
        aligner = OINDiscreteAligner.__new__(OINDiscreteAligner)
        spec = TEMPLATE_SPECS[geo]
        tv = [np.asarray(spec[i]["pos"], dtype=float) for i in sorted(spec)]
        return set(tuple(p) for p in aligner._brute_force_symmetries(tv))

    def test_brute_force_never_exceeds_the_derived_group(self):
        """Brute force may MISS elements, but must never invent one."""
        for geo in GEOMETRY_VERTICES:
            derived = set(geometry_rotation_group(geo))
            extra = self._brute(geo) - derived
            self.assertFalse(extra, f"{geo}: brute force found non-rotations {extra}")

    def test_pbp_is_the_one_disagreement(self):
        disagree = {
            geo
            for geo in GEOMETRY_VERTICES
            if self._brute(geo) != set(geometry_rotation_group(geo))
        }
        self.assertEqual(disagree, {"PBP"}, "the set of geometries where the two differ moved")
        self.assertEqual(len(self._brute("PBP")), 2)
        self.assertEqual(len(geometry_rotation_group("PBP")), 10)


# Colored-vertex helpers for the signature tests. OCT vertex order is
# 0=+z 1=-z 2=+x 3=-x 4=+y 5=-y, so {0,1} {2,3} {4,5} are the three trans pairs.
_N = ("ppy_body", "N", "")
_C = ("ppy_body", "C", "")
FAC = {0: _N, 2: _N, 4: _N, 1: _C, 3: _C, 5: _C}  # three N mutually cis -> facial
MER = {0: _N, 1: _N, 2: _N, 3: _C, 4: _C, 5: _C}  # 0/1 are trans -> meridional


class TestLexMinSignature(unittest.TestCase):
    def test_perm_actually_achieves_the_signature(self):
        for geo, vcolor in (("OCT", FAC), ("OCT", MER), ("SPL", {0: _N, 2: _C})):
            sig, perm = lexmin_vertex_signature(geo, vcolor)
            rebuilt = [("~", "~", "~")] * len(sig)
            for slot, color in vcolor.items():
                rebuilt[perm[slot] if slot < len(perm) else slot] = color
            self.assertEqual(tuple(rebuilt), sig, f"{geo}: perm does not produce the signature")

    def test_fac_and_mer_do_not_collapse(self):
        """The over-folding guard. If this fails, the key silently merges two real isomers."""
        self.assertNotEqual(
            lexmin_vertex_signature("OCT", FAC)[0],
            lexmin_vertex_signature("OCT", MER)[0],
        )

    def test_rotated_relabeling_of_one_isomer_collapses(self):
        """Conformer drift -- the same isomer relabeled by a proper rotation -- must fold."""
        base = lexmin_vertex_signature("OCT", FAC)[0]
        for perm in geometry_rotation_group("OCT"):
            rotated = {perm[slot]: color for slot, color in FAC.items()}
            self.assertEqual(
                lexmin_vertex_signature("OCT", rotated)[0],
                base,
                f"perm {perm} changed the signature",
            )

    def test_tie_break_is_deterministic_and_lex_smallest(self):
        """Repeat calls agree, and the chosen perm is the lex-smallest achieving one.

        The Y2 wave shipped a tie-break that depended on a stereochemical SIGN, which made
        a token reflection-invariant and silently destroyed the chirality it encoded. The
        tie-break here must be a property of the permutation only.
        """
        sig, perm = lexmin_vertex_signature("OCT", FAC)
        for _ in range(5):
            self.assertEqual(lexmin_vertex_signature("OCT", FAC), (sig, perm))
        achieving = []
        for cand in geometry_rotation_group("OCT"):
            arr = [("~", "~", "~")] * 6
            for slot, color in FAC.items():
                arr[cand[slot]] = color
            if tuple(arr) == sig:
                achieving.append(cand)
        self.assertIn(perm, achieving)
        self.assertEqual(perm, min(achieving))


class TestCanonicalSlotPermutation(unittest.TestCase):
    def test_covers_exactly_the_supplied_slots(self):
        mapping = canonical_slot_permutation("OCT", FAC)
        self.assertEqual(set(mapping), set(FAC))

    def test_is_injective(self):
        mapping = canonical_slot_permutation("OCT", FAC)
        self.assertEqual(len(set(mapping.values())), len(mapping))

    def test_slots_beyond_the_template_map_to_themselves(self):
        """A real structure can carry a coordination number the template does not model."""
        vcolor = dict(FAC)
        vcolor[9] = ("extra", "O", "")
        mapping = canonical_slot_permutation("OCT", vcolor)
        self.assertEqual(mapping[9], 9)

    def test_unknown_geometry_is_identity(self):
        vcolor = {0: _N, 1: _C, 2: _N}
        self.assertEqual(canonical_slot_permutation("NOPE", vcolor), {0: 0, 1: 1, 2: 2})

    def test_derive_rotation_group_is_pure(self):
        """Called directly on arbitrary vertices (Lane 5/6 may do this), no cache surprises."""
        a = derive_rotation_group(GEOMETRY_VERTICES["OCT"])
        b = derive_rotation_group(GEOMETRY_VERTICES["OCT"])
        self.assertEqual(a, b)
        self.assertEqual(len(a), 24)


if __name__ == "__main__":
    unittest.main()
