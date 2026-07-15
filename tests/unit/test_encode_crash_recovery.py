"""Guard tests for the two XYZ->OIN encoder crash-recovery paths (S3, v0.4.2).

Both defects made the *input* encode raise before any string was produced -- the
whole molecule was unrepresentable, not merely mis-encoded. Each test drives the
real ``XYZToSMILES.convert`` on a committed fixture distilled from the tmCAT/tmPHOTO
sweep and asserts a clean, stable encode. Both fail against pre-fix code, which
raises a ``ValueError``.

* ``IROXET`` (encode_crash_other): a large saturated ligand whose extended-Huckel
  charge proposal was several electrons off, steering ``_select_lig_mol``'s +-2
  ladder away from any perceivable charge so ``get_lig_mol`` returned ``None``.
  ``get_lig_mol`` now falls through to the wider ``_rescue_unusable_perception``
  sweep, which reaches a usable charge.
* ``JOTJEK`` (kekulize_encode_crash): a fused/charged ring system on which
  ``AC2mol`` left stale aromatic flags, so the first ``SanitizeMol`` could not
  kekulize even though no ring is genuinely quinoid (``stuck_ring_atoms`` empty).
  ``kekulize_safe_sanitize`` now retries with a fresh aromaticity re-perception
  before reporting an honest limitation.

Assertions are kept rdkit-version robust (no exact ``/``\\`` bond-direction
strings): a non-empty encode led by the correct metal token, plus forward-encode
stability (encode twice -> identical), which is the property the round trip needs.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import RDLogger

from oinsmiles.core.translator import XYZToSMILES

RDLogger.DisableLog("rdApp.*")

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))


class TestWidenedChargeSweep(unittest.TestCase):
    """get_lig_mol must not crash when the Huckel charge proposal is far off."""

    def test_saturated_ligand_encodes(self):
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "IROXET_comp_0.xyz"))
        self.assertTrue(oin)
        self.assertTrue(oin.startswith("[Rh"), oin[:40])

    def test_forward_encode_is_stable(self):
        path = os.path.join(_FIXTURES, "IROXET_comp_0.xyz")
        self.assertEqual(XYZToSMILES().convert(path), XYZToSMILES().convert(path))


class TestKekulizeReperception(unittest.TestCase):
    """A stale-aromatic-flag ring system must re-perceive instead of crashing."""

    def test_stale_flag_ligand_encodes(self):
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "JOTJEK_comp_0.xyz"))
        self.assertTrue(oin)
        self.assertTrue(oin.startswith("[Pd"), oin[:40])

    def test_forward_encode_is_stable(self):
        path = os.path.join(_FIXTURES, "JOTJEK_comp_0.xyz")
        self.assertEqual(XYZToSMILES().convert(path), XYZToSMILES().convert(path))


if __name__ == "__main__":
    unittest.main()
