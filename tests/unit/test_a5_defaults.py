"""A5 (v0.4.3) promote-to-default guards.

Pins the two defaults A5's four-arm A/B flipped, so a silent revert is caught:

* the whole-complex vdW acceptance term is ON by default (opt out with
  ``OIN_VDW_ACCEPTANCE=0``) -- the A/B measured clash 92.5%->5.1% with zero
  coordination-number regressions on the improved v0.4.3 pool; and
* ``OIN3DGenerator`` defaults to ``optimizer="ff"`` -- FF + the vdW gate beat the
  g-xTB relax on both clash and round-trip, deterministically and fast (g-xTB stays
  available opt-in).

Both assertions fail against the pre-A5 code (gate off / optimizer="xtb").
"""

import importlib
import inspect
import os
import unittest
from unittest import mock

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generator3d import clash


class TestA5Defaults(unittest.TestCase):
    def test_vdw_gate_on_by_default(self):
        try:
            with mock.patch.dict(os.environ, {}, clear=False):
                os.environ.pop("OIN_VDW_ACCEPTANCE", None)
                importlib.reload(clash)
                self.assertTrue(
                    clash.VDW_ACCEPTANCE_ENABLED,
                    "vdW acceptance must be ON by default (A5 promotion)",
                )
            with mock.patch.dict(os.environ, {"OIN_VDW_ACCEPTANCE": "0"}):
                importlib.reload(clash)
                self.assertFalse(
                    clash.VDW_ACCEPTANCE_ENABLED,
                    "OIN_VDW_ACCEPTANCE=0 must opt out of the gate",
                )
        finally:
            importlib.reload(clash)  # restore module-default state

    def test_engine_default_optimizer_is_ff(self):
        default = inspect.signature(OIN3DGenerator.__init__).parameters["optimizer"].default
        self.assertEqual(default, "ff", "OIN3DGenerator must default to the FF path (A5)")


if __name__ == "__main__":
    unittest.main()
