"""Unit tests for the v0.4.4 fac/mer-aware canonical key (SL0).

The pre-v0.4.4 key renumbered slots by first appearance, so fac-Ir(ppy)3 and
mer-Ir(ppy)3 -- three identical ppy fragments differing only in which polyhedron
vertices the C/N donors occupy -- collapsed to the same key. v0.4.4 replaces the
renumber with a symmetry-canonical vertex signature (lexicographically-minimal
coloring over the coordination polyhedron's proper-rotation group). These tests pin:

* fac != mer at BOTH key layers, while benign slot/fragment permutations stay equal;
* the derived proper-rotation group orders (the signature's core machinery);
* that the local ``_GEOMETRY_VERTICES`` table stays in lockstep with the encoder's
  ``oin_aligner.TEMPLATE_SPECS`` (which is where the slot->vertex convention lives).

Metal ``@SPn`` chirality is DEFERRED, not in the key -- see
``spec/handoffs/v0.4.4/geometry-canonical-slot-key.md`` §4.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.oin.compare import (  # noqa: E402
    _GEOMETRY_VERTICES,
    _geometry_rotation_group,
    canonical_roundtrip_key,
    normalize_oin_for_comparison,
    winding_canonical_key,
)

# Committed goldens (byte-for-byte with test_isomer_divergence / test_regression_stability).
FAC_IRPPY3 = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
MER_IRPPY3 = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"
CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
TRANSPLATIN = "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"

# ACEPUT: a benign slot + fragment-order permutation of the SAME compound (CO/Br swap),
# related by a proper octahedral rotation -- must stay EQUAL under both keys.
ACEPUT_A = "[Re_OCT].C{0}#O.[Br]{1}.Clc1ccn{2}c2c1ccc1ccc3c(Cl)ccn{3}c3c21.C{4}#O.C{5}#O"
ACEPUT_B = "[Re_OCT].[Br]{0}.C{1}#O.Clc1ccn{2}c2c1ccc1ccc3c(Cl)ccn{3}c3c21.C{4}#O.C{5}#O"


def _wk(oin):
    return winding_canonical_key(normalize_oin_for_comparison(oin))


class TestFacMerKey(unittest.TestCase):
    def test_fac_mer_diverge_both_keys(self):
        self.assertNotEqual(_wk(FAC_IRPPY3), _wk(MER_IRPPY3))
        self.assertNotEqual(
            canonical_roundtrip_key(FAC_IRPPY3), canonical_roundtrip_key(MER_IRPPY3)
        )

    def test_cis_trans_diverge_both_keys(self):
        self.assertNotEqual(_wk(CISPLATIN), _wk(TRANSPLATIN))
        self.assertNotEqual(
            canonical_roundtrip_key(CISPLATIN), canonical_roundtrip_key(TRANSPLATIN)
        )

    def test_benign_slot_permutation_stays_equal(self):
        """A rotation-related CO/Br slot + fragment-order swap is the same compound."""
        self.assertEqual(_wk(ACEPUT_A), _wk(ACEPUT_B))
        self.assertEqual(canonical_roundtrip_key(ACEPUT_A), canonical_roundtrip_key(ACEPUT_B))

    def test_keys_are_hashable(self):
        """The v0.4.4 keys are tuples-of-tuples (hashable), unlike the old (str, list)."""
        for key in (_wk(FAC_IRPPY3), canonical_roundtrip_key(FAC_IRPPY3)):
            self.assertIsInstance(hash(key), int)


class TestRotationGroups(unittest.TestCase):
    def test_group_orders(self):
        expected = {
            "LIN": 2,
            "TPL": 6,
            "SPL": 8,
            "TET": 12,
            "TBP": 6,
            "SPY": 4,
            "OCT": 24,
            "PBP": 10,
        }
        for geo, order in expected.items():
            self.assertEqual(
                len(_geometry_rotation_group(geo)), order, f"{geo} rotation-group order"
            )

    def test_permutations_are_valid_bijections(self):
        for geo, verts in _GEOMETRY_VERTICES.items():
            group = _geometry_rotation_group(geo)
            n = len(verts)
            for perm in group:
                self.assertEqual(sorted(perm), list(range(n)), f"{geo} perm not a bijection")

    def test_unknown_geometry_is_none(self):
        self.assertIsNone(_geometry_rotation_group("ZZZ"))


class TestVertexTableInLockstep(unittest.TestCase):
    def test_local_table_matches_template_specs(self):
        """``compare._GEOMETRY_VERTICES`` must track ``oin_aligner.TEMPLATE_SPECS`` directions.

        Guards against silent drift of the slot->vertex convention (the duplication is
        deliberate, to keep compare.py's import graph light -- see its module docstring).
        Compare unit-direction vectors per slot, tolerant to the TET column-normalization
        the aligner applies.
        """
        from oinsmiles.utils.oin_aligner import TEMPLATE_SPECS

        for geo, verts in _GEOMETRY_VERTICES.items():
            self.assertIn(geo, TEMPLATE_SPECS, f"{geo} missing from TEMPLATE_SPECS")
            spec = TEMPLATE_SPECS[geo]
            self.assertEqual(len(verts), len(spec), f"{geo} vertex count")
            for slot, pos in enumerate(verts):
                local = np.asarray(pos, float)
                ref = np.asarray(spec[slot]["pos"], float)
                local = local / np.linalg.norm(local)
                ref = ref / np.linalg.norm(ref)
                self.assertTrue(
                    np.allclose(local, ref, atol=1e-4),
                    f"{geo} slot {slot}: {local} != {ref}",
                )


if __name__ == "__main__":
    unittest.main()
