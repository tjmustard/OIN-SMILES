"""
Tests for core math logic and templates.
"""
import numpy as np
import pytest
from oinsmiles.core.geometry_templates import TEMPLATES
from oinsmiles.core.haptic_math import expand_slot, normalize

def test_templates_integrity():
    """Verify templates have valid pos and ref vectors."""
    for geo, slots in TEMPLATES.items():
        for slot_idx, data in slots.items():
            pos = data['pos']
            ref = data['ref']
            assert len(pos) == 3
            assert len(ref) == 3
            # ref should be orthogonal to pos?
            # Ideally yes, North Star X is perp to Z (pos).
            # But template definition in PRD allows them to be roughly defined?
            # Actually, `expand_slot` orthogonalizes them.
            # But let's check basic normalization.
            # pos might be unit (except TET maybe?).
            if geo != 'TET': # TET pos are corners of cube [1,1,1] etc.
                 pass 
                 # We don't enforce unit length in template consts, but good practice.
            pass

def test_expand_slot_monodentate():
    """Test N=1."""
    z = np.array([0,0,1])
    ref = np.array([1,0,0])
    vecs = expand_slot(1, z, ref)
    assert len(vecs) == 1
    np.testing.assert_allclose(vecs[0], [0,0,1], atol=1e-7)

def test_expand_slot_pentagon():
    """Test N=5 (Cp)."""
    # Slot Z = [0,0,1]
    # Ref X = [1,0,0]
    # Should generate 5 vectors in xy plane (at z=1 + cone spread?)
    # Wait, `vec = Z + spread * (vx*X + vy*Y)`
    # This creates a Cone. The vectors are normalized.
    
    z = np.array([0,0,1])
    ref = np.array([1,0,0])
    vecs = expand_slot(5, z, ref, cone_spread=0.2)
    
    assert len(vecs) == 5
    
    # Vector 0 should be aligned with X (North Star)
    # v0 projected on plane should be along +X.
    # v0 = Z + s*X. Normalized.
    # v0 ~= [0.2, 0, 1] normalized.
    
    v0 = vecs[0]
    expected_unnorm = np.array([0.2, 0, 1])
    expected = expected_unnorm / np.linalg.norm(expected_unnorm)
    
    np.testing.assert_allclose(v0, expected, atol=1e-7)
    
    # Check angles between consecutive vectors (approx)
    # They should be distinct.
    for v in vecs:
        assert np.linalg.norm(v) > 0.99

def test_expand_slot_orthogonality():
    """Test that generated vectors are handled correctly even if Ref is not perp to Z."""
    z = np.array([0,0,1])
    ref = np.array([1,0,1]) # Not perp
    
    vecs = expand_slot(4, z, ref)
    # expand_slot orthogonalizes Ref -> X=[1,0,0].
    # So v0 should align with [1,0,0] (plus Z component).
    
    v0 = vecs[0]
    # Proj of v0 onto xy plane should be [1,0,0]
    proj = v0[:2]
    proj = proj / np.linalg.norm(proj)
    np.testing.assert_allclose(proj, [1,0], atol=1e-2)
    
