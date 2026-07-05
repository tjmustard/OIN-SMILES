"""Unit pin for TASK-41 (WS-1): get_tmc_mol raises a descriptive ValueError
instead of returning a bare None on ligand-perception failure.

Root cause E3 (see spec/worklog/ROUNDTRIP-eta-recovery-handoff.md): when
get_lig_mol cannot build a ligand fragment, get_tmc_mol used to `return None`,
and the sole convert-path caller (core/translator.py) unpacks a 2-tuple, so the
failure surfaced as the opaque "cannot unpack non-iterable NoneType object"
(TiCat3/4). The fixture is a captured TiCat3 *generated* structure whose
clashing geometry over-connects a ligand fragment so perception fails.
"""

import unittest
from pathlib import Path

from oinsmiles.utils.xyz2mol import get_tmc_mol

_FIXTURE = Path(__file__).parent.parent / "fixtures" / "ticat3_generated_broken.xyz"


class TestXyz2MolErrors(unittest.TestCase):
    @unittest.skipUnless(_FIXTURE.exists(), f"fixture missing: {_FIXTURE}")
    def test_get_tmc_mol_raises_valueerror_on_unbuildable_ligand(self):
        """A ligand fragment that cannot be perceived raises a descriptive
        ValueError (not a bare-None-induced TypeError)."""
        with self.assertRaises(ValueError) as ctx:
            get_tmc_mol(_FIXTURE, 0, with_stereo=False)
        self.assertIn("get_lig_mol failed", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
