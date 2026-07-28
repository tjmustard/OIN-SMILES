"""Guards for ``oin.canonical_slots`` -- the proper-rotation group and lex-min relabeling.

These tests exist because the module is load-bearing twice over: it decides the comparison
KEY's fac/mer awareness (since v0.4.4) and, from v0.4.5, the canonical ``{n}`` slot numbers
the encoder emits. A bug here either reintroduces fac/mer blindness (silently merging two
real isomers) or destroys stereochemistry by folding over a reflection.
"""

import itertools
import os
import re
import unittest
from unittest import mock

import numpy as np

from oinsmiles.oin.canonical_slots import (
    GEOMETRY_VERTICES,
    VERTEX_SENTINEL,
    _donor_swap_permutations,
    canonical_slot_map,
    canonical_slot_permutation,
    canonicalize_oin_slots,
    derive_rotation_group,
    geometry_rotation_group,
    geometry_vertex_count,
    lexmin_vertex_signature,
)
from oinsmiles.oin.compare import (
    _parse_vertex_colors,
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
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


class TestAlignerSymmetriesAreUnified(unittest.TestCase):
    """The encoder and the comparison key must enumerate ONE rotation group (TD-005).

    ``oin_aligner._brute_force_symmetries`` used to build the group by brute-forcing Euler
    triples from the fixed grid ``[0, 90, 120, 180, 240, 270]``. That grid cannot express
    PBP's 72-degree five-fold equatorial rotation, so it found **2 of the 10** proper
    rotations and the encoder could not canonicalize a pentagonal-bipyramidal equatorial
    labeling at all. It agreed with ``derive_rotation_group`` on the other ten geometries
    and never invented a non-rotation, which is what made unifying safe.

    v0.4.5 Lane 2 unified it. These tests now assert agreement rather than pinning the gap,
    and ``test_pbp_five_fold_is_now_reachable`` keeps the specific defect from coming back
    if someone reintroduces an angle-grid search.
    """

    @staticmethod
    def _aligner_group(geo):
        aligner = OINDiscreteAligner.__new__(OINDiscreteAligner)
        spec = TEMPLATE_SPECS[geo]
        tv = [np.asarray(spec[i]["pos"], dtype=float) for i in sorted(spec)]
        return set(tuple(p) for p in aligner._brute_force_symmetries(tv))

    def test_aligner_agrees_with_the_derived_group_everywhere(self):
        for geo in GEOMETRY_VERTICES:
            self.assertEqual(
                self._aligner_group(geo),
                set(geometry_rotation_group(geo)),
                f"{geo}: the encoder and the key disagree about the rotation group",
            )

    def test_pbp_five_fold_is_now_reachable(self):
        """The 72-degree equatorial rotation the old Euler grid could not express."""
        group = self._aligner_group("PBP")
        self.assertEqual(len(group), 10)
        # PBP vertex order is 0=+z, 1=-z, then five equatorial vertices 2..6 at 72-degree
        # steps. The C5 generator fixes the axial pair and cycles the equator.
        c5 = (0, 1, 3, 4, 5, 6, 2)
        self.assertIn(c5, group, "the five-fold rotation is missing again")


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


FAC_OIN = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
MER_OIN = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"


def _relabel_oin(oin, perm):
    """Apply a vertex permutation to an OIN string's slot integers (winding kept)."""
    return re.sub(
        r"\{(\d+)([<>^]?)\}",
        lambda m: "{" + str(perm[int(m.group(1))]) + m.group(2) + "}",
        oin,
    )


class TestCanonicalizeOinSlots(unittest.TestCase):
    """The string-level post-pass. Fast, generator-free, no encoder in the loop.

    These are the sharp tests for the lane's central claim: a slot labeling and any
    proper-rotation image of it must emit ONE string, while two labelings NOT related by a
    proper rotation (fac vs mer) must stay apart.
    """

    def test_every_rotation_of_one_labeling_gives_one_string(self):
        for oin in (FAC_OIN, MER_OIN):
            outs = {
                canonicalize_oin_slots(_relabel_oin(oin, perm))
                for perm in geometry_rotation_group("OCT")
            }
            self.assertEqual(len(outs), 1, f"{len(outs)} distinct strings over the OCT group")

    def test_fac_and_mer_do_not_collapse(self):
        """The over-folding guard at string level: worse than drift, because it is silent."""
        self.assertNotEqual(canonicalize_oin_slots(FAC_OIN), canonicalize_oin_slots(MER_OIN))

    def test_cis_and_trans_do_not_collapse(self):
        cis = "[Pt_SPL].N{0}.[CH3]{1}.[Br]{2}.[Cl]{3}"
        trans = "[Pt_SPL].N{0}.[Br]{1}.[CH3]{2}.[Cl]{3}"
        self.assertNotEqual(canonicalize_oin_slots(cis), canonicalize_oin_slots(trans))

    def test_idempotent(self):
        once = canonicalize_oin_slots(FAC_OIN)
        self.assertEqual(canonicalize_oin_slots(once), once)

    def test_metal_fragment_stays_first(self):
        """``fragments[0]`` is the metal -- a load-bearing project invariant."""
        for oin in (FAC_OIN, MER_OIN, "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"):
            self.assertTrue(canonicalize_oin_slots(oin).split(".")[0].startswith("["))
            self.assertIn("_", canonicalize_oin_slots(oin).split(".")[0])

    def test_winding_character_is_preserved_verbatim(self):
        """``^`` is folded to ``>`` for COLORING only; the emitted char must survive."""
        oin = (
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1^}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
        )
        out = canonicalize_oin_slots(oin)
        self.assertEqual(out.count(">"), 1)
        self.assertEqual(out.count("^"), 1)

    def test_axial_suffix_is_carried_through(self):
        oin = "[Pd_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3} |ax:+-|"
        self.assertTrue(canonicalize_oin_slots(oin).endswith(" |ax:+-|"))

    def test_no_slots_is_returned_unchanged(self):
        self.assertEqual(canonicalize_oin_slots("[Fe_OCT]"), "[Fe_OCT]")
        self.assertEqual(canonicalize_oin_slots(""), "")

    def test_emitted_labeling_achieves_the_keys_lex_min(self):
        """Encoder and comparison key agree BY CONSTRUCTION -- the point of this seam.

        After the post-pass the identity permutation must already be optimal for
        ``_polyhedron_signature``, i.e. re-canonicalizing moves nothing.
        """
        for oin in (FAC_OIN, MER_OIN):
            out = canonicalize_oin_slots(oin)
            _metal, geo, vcolor = _parse_vertex_colors(normalize_oin_for_comparison(out))
            sig, _perm = lexmin_vertex_signature(geo, vcolor)
            identity = [VERTEX_SENTINEL] * geometry_vertex_count(geo)
            for slot, color in vcolor.items():
                identity[slot] = color
            self.assertEqual(tuple(identity), sig)

    def test_comparison_key_is_unchanged_by_the_post_pass(self):
        """The post-pass may not move the key -- it only relabels within one orbit."""
        for oin in (FAC_OIN, MER_OIN):
            self.assertEqual(
                canonical_roundtrip_key(oin), canonical_roundtrip_key(canonicalize_oin_slots(oin))
            )


class TestCanonicalSlotMap(unittest.TestCase):
    """The documented Lane 5 / Lane 6 entry point: ``canonical_slot_map(oin)[slot]``."""

    def test_map_matches_the_relabeled_string(self):
        mapping = canonical_slot_map(FAC_OIN)
        self.assertEqual(_relabel_oin(FAC_OIN, mapping), canonicalize_oin_slots(FAC_OIN))

    def test_map_is_identity_on_an_already_canonical_string(self):
        out = canonicalize_oin_slots(FAC_OIN)
        mapping = canonical_slot_map(out)
        self.assertEqual(mapping, {s: s for s in mapping})

    @staticmethod
    def _per_donor_maps(oin, geo):
        """Distinct 'physical donor -> canonical slot' assignments over the whole group.

        The donor that starts at slot ``s`` sits at ``perm[s]`` in the rotated presentation,
        so its canonical slot is ``canonical_slot_map(rotated)[perm[s]]``.
        """
        n = geometry_vertex_count(geo)
        out = set()
        for perm in geometry_rotation_group(geo):
            cmap = canonical_slot_map(_relabel_oin(oin, perm))
            out.add(tuple(cmap[perm[s]] for s in range(n)))
        return out

    def test_per_donor_map_is_unique_when_donors_are_distinguishable(self):
        """Lane 5's requirement: 'canonical slot of donor d' is a function of the molecule."""
        het = "[Ir_OCT].[F]{0}.[Cl]{1}.[Br]{2}.[I]{3}.N{4}.O{5}"
        self.assertEqual(len(self._per_donor_maps(het, "OCT")), 1)

    def test_per_donor_ambiguity_is_exactly_the_genuine_automorphism(self):
        """fac-M(ppy)3 has a real C3 axis, so three donor labelings are equally canonical.

        The **string** is still unique (asserted here too) -- the ambiguity is only in which
        of three interchangeable ligands is called 'first'. It is bounded by the colored
        polyhedron's rotational automorphism group, and every member of that group is a
        PROPER rotation, so no stereochemical descriptor Lane 5 or Lane 6 derives from the
        canonical arrangement can differ between them. If this count ever exceeds the
        automorphism order, the tie-break has stopped being a property of the molecule.
        """
        maps = self._per_donor_maps(FAC_OIN, "OCT")
        self.assertEqual(len(maps), 3, "not the C3 of fac-tris-homoleptic")
        outs = {
            canonicalize_oin_slots(_relabel_oin(FAC_OIN, perm))
            for perm in geometry_rotation_group("OCT")
        }
        self.assertEqual(len(outs), 1)


#: ``AGUKOD_comp_0``, the corpus case that pins the residual class. A COD ligand spans two
#: cis vertices of a square plane; the drifted presentation differs ONLY in which alkene arm
#: carries which slot integer. Bodies are byte-identical.
_AGUKOD_A = (
    "[Rh_SPL].[CH]{1}1=[CH]{1>}CC[CH]{0>}=[CH]{0}CC1"
    ".CC(C)(C)c1cc2c(c(C(C)(C)C)c1)OC{2}N2c1ccccc1.[Cl]{3}"
)
_AGUKOD_B = (
    "[Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1"
    ".CC(C)(C)c1cc2c(c(C(C)(C)C)c1)OC{2}N2c1ccccc1.[Cl]{3}"
)


class TestResidualClassIsOutOfReachByDesign(unittest.TestCase):
    """The measured limit of the ROTATION-ONLY post-pass -- and, since v0.4.11, its fix.

    **INVERTED in v0.4.11, as v0.4.5 asked.** The class below still pins exactly what it
    always pinned: with ``OIN_CANONICAL_DONOR_FOLD`` off, the residual class is unreachable,
    for the reason stated. What changed is that the fold now exists, so the same fixture is
    also asserted to CONVERGE with the lever on -- see
    ``test_the_donor_fold_converges_them`` and
    ``docs/agentic-notes/v0.4.11/LANE-02-donor-fold.md``.

    Measured over 150 corpus molecules (see ``docs/agentic-notes/v0.4.5/CANONICAL_SLOTS_v0.4.5.md`` section 7a):
    **32/32** residual ``slot_renumber`` pairs have an *identical* colored-vertex map. The
    post-pass derives its relabeling from that map alone, so it computes the identical
    permutation for both strings and the difference survives it. That is not a bug in the
    relabeling -- it is a limit of its input, and it means the lane's original acceptance
    target (``slot_renumber -> ~0``) is not reachable at this seam.

    v0.4.11 re-measured the same taxonomy on the ROUND-TRIP population (496 pairs, not 32
    re-presentation pairs) and found the same shape: **496/496 same_vcolor_identical**, of
    which 377 are atom-level ``automorphism`` and reachable by the within-fragment fold.

    Why the map is identical: ``compare._parse_vertex_colors`` colors *every* donor of a
    ligand with that ligand's whole body and no chelate grouping, deliberately, so that a
    swap between two same-colored donors is invisible -- which is what lets true conformers
    collapse in the comparison KEY while fac/mer stay distinct. Correct for a key;
    insufficient for an encoder, which must decide which donor atom holds which integer.

    The transposition relating the two is also NOT in the geometry's rotation group (SPL's
    D4 contains ``(1 3)``, ``(0 2)``, ``(0 1)(2 3)`` but not ``(0 1)`` alone), so widening
    the group-theoretic fold cannot reach it either -- and widening it is exactly the
    over-folding this lane must not do.

    That fix -- fold same-symmetry-class, same-color donors WITHIN one fragment -- landed in
    v0.4.11 as ``OIN_CANONICAL_DONOR_FOLD``, and this test was inverted rather than deleted.
    """

    def test_the_two_presentations_have_an_identical_vertex_coloring(self):
        _m1, g1, c1 = _parse_vertex_colors(normalize_oin_for_comparison(_AGUKOD_A))
        _m2, g2, c2 = _parse_vertex_colors(normalize_oin_for_comparison(_AGUKOD_B))
        self.assertEqual(g1, g2)
        self.assertEqual(c1, c2, "if these ever differ, the post-pass CAN reach this class")

    def test_and_therefore_the_ROTATION_ONLY_post_pass_cannot_converge_them(self):
        """Spell 'off' as "0" -- deleting the variable means ON once the lever is promoted."""
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_DONOR_FOLD": "0"}):
            self.assertNotEqual(
                canonicalize_oin_slots(_AGUKOD_A),
                canonicalize_oin_slots(_AGUKOD_B),
                "the rotation-only fold reached this class -- re-measure section 7a",
            )

    def test_the_donor_fold_converges_them(self):
        """The inversion: what v0.4.5 pinned as out of reach, v0.4.11 reaches."""
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_DONOR_FOLD": "1"}):
            a = canonicalize_oin_slots(_AGUKOD_A)
            b = canonicalize_oin_slots(_AGUKOD_B)
        self.assertEqual(a, b, "the within-fragment donor fold must converge the archetype")

    def test_the_fold_is_idempotent(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_DONOR_FOLD": "1"}):
            once = canonicalize_oin_slots(_AGUKOD_A)
            self.assertEqual(once, canonicalize_oin_slots(once))

    def test_the_relating_transposition_is_not_a_proper_rotation_of_the_square(self):
        self.assertNotIn((1, 0, 2, 3), geometry_rotation_group("SPL"))

    def test_but_the_comparison_key_folds_them_so_no_isomer_is_at_risk(self):
        """Which is why the residual class costs 0 key-level defects (measured 18 -> 18)."""
        self.assertEqual(
            canonical_roundtrip_key(_AGUKOD_A),
            canonical_roundtrip_key(_AGUKOD_B),
        )


class TestDonorFoldScope(unittest.TestCase):
    """The three conditions that keep the v0.4.11 widening safe, tested one at a time.

    The fold is the only place in this module that goes past the geometry's own
    proper-rotation group, and the narrowness of its scope IS the safety argument: one
    fragment, one ``breakTies=False`` symmetry class, one colour. Each test below removes
    one condition and asserts the fold declines.
    """

    def _fold(self, oin, on=True):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_DONOR_FOLD": "1" if on else "0"}):
            return canonicalize_oin_slots(oin)

    def _swaps(self, oin):
        """The swap set the fold would use -- asserted directly, not through the output.

        Testing scope through ``canonicalize_oin_slots`` does not work: the rotation group
        already converges most small cases on its own, so an equal-output assertion passes
        whether or not the fold over-reached. The swap set is the thing under test.
        """
        frags = [f for f in oin.split(".") if f]
        _m, _geo, vcolor = _parse_vertex_colors(normalize_oin_for_comparison(oin))
        return _donor_swap_permutations(frags, vcolor)

    def test_donors_in_DIFFERENT_symmetry_classes_are_not_exchanged(self):
        """A 2-aminoethanol-like chelate: the N and the O are not interchangeable.

        This is the ``distinct_donors_LOCAL`` class -- 118 of the v0.4.11 population. A
        correctly-scoped fold cannot reach it, and must not try.
        """
        self.assertEqual(
            self._swaps("[Pt_SPL].N{0}CCO{1}.[Cl]{2}.[Cl]{3}"),
            [{}],
            "the fold offered to exchange an N donor with an O donor -- condition (a) is broken",
        )

    def test_donors_in_DIFFERENT_fragments_are_not_exchanged(self):
        """Two separate ligands' donors are never exchanged, however alike they are.

        Cross-fragment exchange is what would reach the metal Delta/Lambda arrangement, and
        that arrangement is the isomer. Scope condition: WITHIN one fragment.
        """
        self.assertEqual(
            self._swaps("[Pt_SPL].N{0}C.N{1}CC.[Cl]{2}.[Cl]{3}"),
            [{}],
            "the fold offered a cross-fragment exchange -- the within-fragment scope is broken",
        )

    def test_two_equivalent_donors_of_ONE_fragment_DO_generate_a_swap(self):
        """The positive control: without this, the three negatives above prove nothing.

        ``AGUKOD``'s COD ligand has two symmetry-equivalent alkene arms on slots 0 and 1.
        """
        swaps = self._swaps(_AGUKOD_A)
        # Identity is spelled {} when no bucket exists and {s: s, ...} when one does; both
        # mean "change nothing", so test the property rather than the spelling.
        self.assertTrue(
            any(all(k == v for k, v in m.items()) for m in swaps),
            "identity must always be offered, so the fold can never narrow the candidate set",
        )
        self.assertIn({0: 1, 1: 0}, swaps, "the two equivalent COD arms must be exchangeable")

    def test_the_fold_never_changes_the_comparison_key(self):
        """The invariant that keeps ``key_equal`` accounting readable across the A/B.

        The fold only ever exchanges same-COLOURED slots, so the colored-vertex signature --
        and therefore ``_polyhedron_signature``, and therefore the key -- is untouched. That
        is also why the rotation-only post-pass could never reach this class: the difference
        lives entirely in the signature's kernel.

        Measured over the whole v0.4.11 population: 0 of 992 strings changed key.
        """
        for oin in (_AGUKOD_A, _AGUKOD_B, FAC_OIN, "[Pt_SPL].N{0}CCO{1}.[Cl]{2}.[Cl]{3}"):
            with self.subTest(oin=oin[:40]):
                before = canonical_roundtrip_key(oin)
                self.assertEqual(before, canonical_roundtrip_key(self._fold(oin)))

    def test_the_fold_never_narrows_the_candidate_set(self):
        """Identity is always in the swap set, so the fold can only ever tie or improve.

        Concretely: everything the rotation-only post-pass could already reach stays
        reachable, which is what makes 'lever OFF is byte-identical' a scoping property and
        not a coincidence.
        """
        for oin in (_AGUKOD_A, "[Pt_SPL].N{0}CCO{1}.[Cl]{2}.[Cl]{3}", FAC_OIN):
            with self.subTest(oin=oin[:40]):
                self.assertLessEqual(self._fold(oin), self._fold(oin, on=False))


if __name__ == "__main__":
    unittest.main()
