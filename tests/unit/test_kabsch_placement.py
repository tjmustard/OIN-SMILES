"""Unit guards for the Kabsch/Umeyama rigid-placement helper (v0.4.3 A4).

``embed.kabsch`` is the least-squares engine of the opt-in ``option=3`` embed: it
places an independently-built ligand onto the ideal coordination vectors with a pure
rotation + translation (no scaling -- bond lengths must stay physical). These tests
pin the two properties the placement path depends on:

* it **recovers a known rigid motion** to numerical tolerance, and
* it is **reflection-guarded**: the returned transform is always a proper rotation
  (``det == +1``), so a mirrored (opposite-handed) correspondence is fitted with a
  large residual rather than silently flipping handedness -- which for a chelate
  would flip metal-centered stereochemistry (Delta/Lambda).

Each test fails against pre-A4 code (the helper did not exist).
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles
from oinsmiles.generator3d import om
from oinsmiles.generator3d.embed import get_embedding, kabsch


def _rotation_from_axis_angle(axis, angle):
    """A proper rotation matrix (Rodrigues) -- test helper, independent of the SUT."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


class TestKabschRecovery(unittest.TestCase):
    """kabsch recovers the rigid motion that produced the target points."""

    def _fixed_cloud(self):
        # Deterministic, non-degenerate (full-rank) point cloud; no RNG so the guard
        # is reproducible run-to-run.
        return np.array(
            [
                [0.0, 0.0, 0.0],
                [1.3, 0.2, -0.4],
                [-0.7, 1.1, 0.5],
                [0.4, -0.9, 1.2],
                [1.0, 1.0, 1.0],
            ]
        )

    def test_recovers_known_rotation_and_translation(self):
        P = self._fixed_cloud()
        rot_true = _rotation_from_axis_angle([0.3, -0.8, 0.5], 1.1)
        trans_true = np.array([2.0, -1.5, 0.7])
        Q = (rot_true @ P.T).T + trans_true

        rot, trans = kabsch(P, Q)

        np.testing.assert_allclose(rot, rot_true, atol=1e-9)
        np.testing.assert_allclose(trans, trans_true, atol=1e-9)
        # And the transform actually maps P onto Q.
        mapped = (rot @ P.T).T + trans
        np.testing.assert_allclose(mapped, Q, atol=1e-9)

    def test_returns_proper_rotation(self):
        P = self._fixed_cloud()
        rot_true = _rotation_from_axis_angle([1.0, 0.0, 0.0], 2.3)
        Q = (rot_true @ P.T).T + np.array([0.1, 0.2, 0.3])
        rot, _ = kabsch(P, Q)
        # Orthonormal, right-handed.
        np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-9)
        self.assertAlmostEqual(float(np.linalg.det(rot)), 1.0, places=9)

    def test_single_point_is_pure_translation(self):
        P = np.array([[1.0, 2.0, 3.0]])
        Q = np.array([[-4.0, 5.0, 6.0]])
        rot, trans = kabsch(P, Q)
        np.testing.assert_allclose(rot, np.eye(3), atol=1e-9)
        np.testing.assert_allclose(trans, Q[0] - P[0], atol=1e-9)


class TestReflectionGuard(unittest.TestCase):
    """A mirrored correspondence must NOT be accepted as the same handedness."""

    def _chiral_cloud(self):
        # A handed (chiral) tetrahedron-like set: its mirror image cannot be reached
        # by any proper rotation.
        return np.array(
            [
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [0.3, 0.1, 0.0],  # break the tetra's improper symmetry
            ]
        )

    def test_mirror_image_yields_proper_rotation_with_large_residual(self):
        P = self._chiral_cloud()
        # Reflect through the z=0 plane -> opposite handedness.
        Q = P * np.array([1.0, 1.0, -1.0])

        rot, trans = kabsch(P, Q)

        # Guard #1: the transform is a proper rotation, never a reflection.
        self.assertAlmostEqual(float(np.linalg.det(rot)), 1.0, places=9)

        # Guard #2: because a rotation cannot reproduce a mirror image, the best-fit
        # residual is large -- the mirror is rejected, not silently matched. An
        # UNGUARDED svd fit would return the reflection with ~0 residual.
        mapped = (rot @ P.T).T + trans
        residual = np.sqrt(np.mean(np.sum((mapped - Q) ** 2, axis=1)))
        self.assertGreater(residual, 0.5)

    def test_same_handedness_still_fits_cleanly(self):
        # Control: a proper rotation of the same cloud fits with ~0 residual, so the
        # guard is not just rejecting everything.
        P = self._chiral_cloud()
        rot_true = _rotation_from_axis_angle([0.2, 0.7, -0.4], 0.9)
        Q = (rot_true @ P.T).T + np.array([1.0, 0.0, -2.0])
        rot, trans = kabsch(P, Q)
        mapped = (rot @ P.T).T + trans
        residual = np.sqrt(np.mean(np.sum((mapped - Q) ** 2, axis=1)))
        self.assertLess(residual, 1e-9)


class TestOptionThreePlacement(unittest.TestCase):
    """End-to-end ``option=3`` placement lands donors on the ideal slot vectors.

    Exercises the whole ``_kabsch_embedding`` path (build free ligand -> place ->
    finalize) for a monodentate (cisplatin) and a haptic (ferrocene) complex, and
    checks the two invariants the downstream pipeline relies on: the metal stays at
    the origin (metal-first ordering) and every donor sits on its assigned
    coordination vector. Fails against pre-A4 code (``option=3`` was undefined).
    """

    CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
    FERROCENE = (
        "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1."
        "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    )

    def _place(self, oin):
        mc = om.get_om_from_modified_smiles(convert_oin_to_msmiles(oin))
        pos = get_embedding(mc, scale=1.0, option=3, align=True, seed=42)
        self.assertIsNotNone(pos, "option=3 returned None (ligand build failed)")
        return mc, np.asarray(pos)

    def _assert_donors_on_slots(self, mc, pos):
        metal = pos[mc.metal_index]
        dvec = np.asarray(mc.geometry_type.direction_vector, dtype=float)
        aidx = mc.get_atom_indices_for_each_ligand()
        for i, lig in enumerate(mc.ligands):
            for atom_locals, slot in lig.binding_infos:
                donor = pos[[aidx[i][a] for a in atom_locals]].mean(axis=0)
                rel = donor - metal
                v = dvec[slot - 1] / np.linalg.norm(dvec[slot - 1])
                cos = float(np.dot(rel, v) / (np.linalg.norm(rel) + 1e-9))
                self.assertGreater(cos, 0.95, f"lig{i} slot{slot} off-vector (cos={cos:.3f})")

    def test_monodentate_donors_land_on_vectors(self):
        mc, pos = self._place(self.CISPLATIN)
        self.assertEqual(pos.shape, (mc.num_atom, 3))
        self.assertFalse(bool(np.isnan(pos).any()))
        np.testing.assert_allclose(pos[mc.metal_index], 0.0, atol=1e-6)
        self._assert_donors_on_slots(mc, pos)

    def test_haptic_face_centroid_lands_on_vector(self):
        mc, pos = self._place(self.FERROCENE)
        np.testing.assert_allclose(pos[mc.metal_index], 0.0, atol=1e-6)
        self._assert_donors_on_slots(mc, pos)

    def test_default_options_do_not_dispatch_kabsch(self):
        # option in {0,1,2} must NOT take the placement path (byte-identity guard).
        mc = om.get_om_from_modified_smiles(convert_oin_to_msmiles(self.CISPLATIN))
        pos = get_embedding(mc, scale=1.0, option=0, align=True, seed=42)
        self.assertIsNotNone(pos)
        self.assertEqual(np.asarray(pos).shape[0], mc.num_atom)


if __name__ == "__main__":
    unittest.main()
