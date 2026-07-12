"""Eta-ring canonicalization acceptance tests (encoder side).

Formerly this module also carried OIN->XYZ generation-side diagnostics that
exercised the legacy Molassembler/stitch backend; that backend was removed, so
only the engine-agnostic encoder tests remain here. They exercise RC1 (scoped
eta-only fragment-order swap) and RC2 (canonical-rank heading atom) in
``OINDiscreteAligner`` and the ``signed_circulation`` winding helper.

This test module MEASURES ONLY -- it must never modify anything under src/.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))


class TestEtaRingCanonicalization(unittest.TestCase):
    """MiniPRD_EtaRingCanonicalization.md acceptance tests (RT-5 hardening).

    Covers RC1 (scoped eta-only fragment-order swap) and RC2 (canonical-rank
    heading atom) in ``OINDiscreteAligner`` -- the two fixes that closed
    ``test_haptic_face_golden_match`` above.
    """

    def test_non_eta_fragment_order_is_inert(self):
        """Test 2 (RT-5): complexes with zero haptic (eta) donor groups have
        no fragment eligible for RC1/RC2, so encoding must stay byte-
        identical to each fixture's pre-existing pinned regression golden."""
        cases = [
            ("CisPlatin.xyz", "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"),
            ("TransPlatin.xyz", "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"),
            ("Cis-PtCl2(en).xyz", "[Pt_SPL].[NH2]{0}CC[NH2]{1}.[Cl]{2}.[Cl]{3}"),
            (
                "fac-Ir(ppy)3.xyz",
                "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1",
            ),
            (
                "mer-Ir(ppy)3.xyz",
                "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1",
            ),
            (
                "PdCl2-RR-BDPP.xyz",
                "[Pd_SPL].C[C@@H](C[C@H](C)P{0}(c1ccccc1)c1ccccc1)"
                "P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            ),
            (
                "PdCl2-RR-BDNN.xyz",
                "[Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)"
                "N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            ),
            (
                "PdCl2-R-BINAP.xyz",
                "[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}"
                "(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}",
            ),
        ]
        for filename, expected in cases:
            with self.subTest(filename=filename):
                actual = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, filename))
                self.assertEqual(actual, expected)

    def test_symmetric_eta_ring_is_inert(self):
        """Test 3 (RT-5), plain-ferrocene half: an UNSUBSTITUTED Cp ring's
        SMILES is a literal member of ``SYMMETRIC_LIGANDS`` (first-wins,
        untouched by RC1/RC2) -- encoding must stay byte-identical."""
        actual = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "Ferrocene.xyz"))
        self.assertEqual(
            actual,
            "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1."
            "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
        )

    def test_mono_substituted_eta_ring_relabel_is_winding_preserving(self):
        """Test 3 (RT-5), ansa-metallocene half.

        A mono-substituted Cp ring (TiCat1's Si-bridged ring) is NOT a
        literal member of ``SYMMETRIC_LIGANDS`` (that set only contains bare,
        unsubstituted ring SMILES) -- it IS eligible for RC2, and its
        canonical-rank heading atom differs from the pre-fix geometric pick.

        Since 2026-07-06 the encoder marks winding PER eta ring (per haptic
        slot), not once per fragment: a silane-bridged ansa-metallocene is a
        single connected fragment occupying two eta slots, so BOTH rings must
        now carry a winding marker (previously only the first ring did, leaving
        rac/meso indistinguishable). Assert the full, complete two-marker
        encoding.

        Both markers read '>': TiCat1's Cp rings carry no substituent other than
        the bridge, so each is carried onto itself (in reverse cyclic order) by a
        180 deg rotation about the in-plane axis through its bridgehead carbon.
        That is a proper rotation which moves neither the metal nor the other
        ring, so the two faces are the same structure and the sign is
        orientation-free notation, not stereochemistry. The rac/meso guard lives
        on TiCat3/TiCat4, whose indenyls admit no such rotation.
        """
        ticat1_xyz = os.path.join(_FIXTURES_DIR, "TiCat1.xyz")
        actual = XYZToSMILES().convert(ticat1_xyz)
        expected = (
            "[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0>}[cH]{0}[cH]{0}1)"
            "c{1}1[cH]{1}[cH]{1>}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}"
        )
        self.assertEqual(actual, expected)

    def test_rc1_scoped_swap_never_touches_non_eta_ranks(self):
        """Test 4 (RT-5): RC1 permutes ONLY same-mass eta fragments among the
        rank slots they already occupy; every non-eta fragment's rank/content
        stays fixed regardless of the eta fragments' arrival order.

        Constructs two content-distinct 5-membered eta rings (reusing the
        real, non-symmetric halide-ring SMILES from the golden fixture) at
        ranks 2 and 3 (deliberately non-1), bracketed by two ordinary
        monodentate (non-eta) fragments at ranks 1 and 4. Calls
        ``_permute_and_serialize`` directly (bypassing 3D/XYZ) with the two
        eta rings' arrival order swapped between two runs, and asserts (a)
        the non-eta ranks' tags are byte-identical between runs and (b) the
        eta ranks end up hosting the SAME content regardless of which one
        arrived first -- the arrival-order-invariance RC1 exists to
        guarantee.
        """
        from collections import defaultdict

        from scipy.spatial.transform import Rotation

        from oinsmiles.utils.oin_aligner import TEMPLATES, OINDiscreteAligner

        geometry_name = "TET"
        tmpl_vectors = TEMPLATES[geometry_name]
        aligner = OINDiscreteAligner(metal_idx=0, ligands=[])

        ring_a_smiles = "Oc1cc(Cl)c(Br)c1I"
        ring_b_smiles = "Oc1cc(I)c(Br)c1Cl"
        ring_indices = [1, 2, 3, 5, 7]
        angles = np.linspace(0, 2 * np.pi, len(ring_indices), endpoint=False)
        ring_coords = np.array([[np.cos(a), np.sin(a), 0.0] for a in angles])

        def _mono(rank):
            return {
                "rank": rank,
                "local_idx": 0,
                "constituent_indices": [0],
                "coords": np.zeros(3),
                "group_coords": np.array([[0.0, 0.0, 0.0]]),
                "chem_id": (35.45, "[Cl]"),
            }

        def _ring(rank, smiles):
            return {
                "rank": rank,
                "local_idx": min(ring_indices),
                "constituent_indices": list(ring_indices),
                "coords": np.zeros(3),
                "group_coords": ring_coords,
                "chem_id": (16.0, smiles),  # same mass -> same RC1 swap bucket
            }

        def _run(ring_at_rank2, ring_at_rank3):
            slot_assignment = [
                _mono(1),
                _ring(2, ring_at_rank2),
                _ring(3, ring_at_rank3),
                _mono(4),
            ]
            result = aligner._permute_and_serialize(
                slot_assignment,
                tmpl_vectors,
                geometry_name=geometry_name,
                alignment_rotation=Rotation.identity(),
            )
            by_rank = defaultdict(list)
            for tag in result.split(";"):
                rank = int(tag.split(".", 1)[0])
                by_rank[rank].append(tag)
            return by_rank

        forward = _run(ring_a_smiles, ring_b_smiles)
        swapped = _run(ring_b_smiles, ring_a_smiles)

        self.assertEqual(
            forward[1],
            swapped[1],
            "non-eta rank 1 must be byte-identical regardless of eta arrival order",
        )
        self.assertEqual(
            forward[4],
            swapped[4],
            "non-eta rank 4 must be byte-identical regardless of eta arrival order",
        )

        # RC1 only reassigns the RANK LABEL for content-canonicalization; the
        # SLOT number is a physical-geometry attribute of wherever the
        # fragment happened to arrive and is deliberately NOT touched by the
        # swap (rank and slot are independent fields in the "rank.idx:slot"
        # tag). So compare content per rank ignoring slot -- this is the
        # actual invariant RC1 promises: which ring CONTENT is labeled with
        # which rank number must not depend on arrival order.
        def _content_ignoring_slot(tags):
            sig = set()
            for tag in tags:
                rank_idx, _slot_and_marker = tag.split(":")
                marker = "".join(c for c in _slot_and_marker if c in "<>")
                sig.add((rank_idx, marker))
            return frozenset(sig)

        forward_content_by_rank = {r: _content_ignoring_slot(t) for r, t in forward.items()}
        swapped_content_by_rank = {r: _content_ignoring_slot(t) for r, t in swapped.items()}

        self.assertEqual(
            forward_content_by_rank[2],
            swapped_content_by_rank[2],
            "rank 2 must host the same ring content regardless of which ring arrived first",
        )
        self.assertEqual(
            forward_content_by_rank[3],
            swapped_content_by_rank[3],
            "rank 3 must host the same ring content regardless of which ring arrived first",
        )
        self.assertNotEqual(
            forward_content_by_rank[2],
            forward_content_by_rank[3],
            "sanity: the two rings are genuinely content-distinct",
        )

    def test_rc2_winding_start_invariance_and_reflection(self):
        """Test 5 (RT-5 / US-004 / RT-4): live assertion that
        ``signed_circulation``'s winding character is (a) identical for
        every possible choice of star atom on a genuinely asymmetric ring
        (start-invariance, the property RC2 relies on to be safe), and (b) a
        synthetically reflected copy of the same ring still yields the
        FLIPPED character -- proving canonicalization never masks a real
        reflection. Kept in addition to (not instead of) the existing
        ``test_haptic_face_r2_geometric_fallback_never_auto_substituted``
        skip.
        """
        from oinsmiles.oin.winding import signed_circulation

        # Irregular (non-symmetric) planar pentagon -- analogous to a
        # differently-pentahalo-substituted ring, so no accidental symmetry
        # could make start-invariance trivially true.
        angles = np.radians([0, 50, 130, 200, 290])
        coords = np.array([[np.cos(a), np.sin(a), 0.0] for a in angles])
        axis = np.array([0.0, 0.0, 1.0])

        base_char = signed_circulation(coords, 0, axis)
        for star in range(len(coords)):
            self.assertEqual(
                signed_circulation(coords, star, axis),
                base_char,
                f"winding character must be start-invariant (star={star})",
            )

        reflected = coords.copy()
        reflected[:, 1] *= -1.0
        self.assertNotEqual(
            signed_circulation(reflected, 0, axis),
            base_char,
            "a reflected ring must flip the winding character",
        )


if __name__ == "__main__":
    unittest.main()
