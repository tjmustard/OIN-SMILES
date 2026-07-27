"""Unit pin for TASK-41 (WS-1): get_tmc_mol raises a descriptive ValueError
instead of returning a bare None on ligand-perception failure.

Root cause E3 (see spec/worklog/ROUNDTRIP-eta-recovery-handoff.md): when
get_lig_mol cannot build a ligand fragment, get_tmc_mol used to `return None`,
and the sole convert-path caller (core/translator.py) unpacks a 2-tuple, so the
failure surfaced as the opaque "cannot unpack non-iterable NoneType object"
(TiCat3/4).

⚠ v0.4.5: THE CONTRACT IS NOW TESTED BY FAULT INJECTION, NOT BY THE FIXTURE.

The original test drove the error path with a captured TiCat3 *generated* structure whose
clashing geometry over-connects a ligand fragment. Promoting ``OIN_STABLE_METAL_AC`` to
default-ON changed that: the lever caps valences highest-Z-first (so that perception no longer
depends on input atom order), which on this degenerate geometry lets the titanium absorb the
contested bonds instead of the walk dead-ending. Perception then *succeeds* -- and returns
nonsense: 48 atoms in 8 fragments, seven bare ``[H+]`` ions and a ``[Ti-14]`` centre.

Two separate facts follow, and both are recorded rather than papered over:

1. The error-path contract still matters, so it is exercised directly below by making
   ``get_lig_mol`` fail. A contract test that depends on a fixture staying unbuildable is one
   perception improvement away from silently testing nothing -- which is exactly what happened
   here, and it would have hidden a regression back to the bare-None TypeError.

2. On a *broken* input, "perceives a nonsense graph" is worse than "fails loudly", so the
   second test pins what the lever actually does today instead of pretending it is fine.
   On real data the lever is clean (capstone A/B: 145 molecules fixed, zero correctness
   regressions; ``geometry_tag_shift`` 0/298 ``[M_XXX]`` changes), so this is a degenerate-input
   concern, not a shipped-accuracy one. The follow-up worth considering is a sanity gate that
   rejects a perceived molecule containing isolated bare-proton fragments; it is NOT added here
   because charged hydrides are legitimate and the gate needs its own corpus A/B first.
"""

import unittest
from pathlib import Path
from unittest import mock

from oinsmiles.utils import perception_tmc as perception_module
from oinsmiles.utils.perception_tmc import get_tmc_mol

_FIXTURE = Path(__file__).parent.parent / "fixtures" / "ticat3_generated_broken.xyz"


class TestXyz2MolErrors(unittest.TestCase):
    @unittest.skipUnless(_FIXTURE.exists(), f"fixture missing: {_FIXTURE}")
    def test_get_tmc_mol_raises_valueerror_when_get_lig_mol_fails(self):
        """An unbuildable ligand fragment raises a descriptive ValueError.

        Fault-injected so the assertion is about ``get_tmc_mol``'s error handling and nothing
        else. The failure must not surface as the bare-None-induced TypeError from
        ``core/translator.py`` unpacking a 2-tuple.
        """
        # (None, charge): get_lig_mol returns a 2-tuple that the call site unpacks before the
        # `if not lig_mol` guard, so injecting a bare None would fail in the unpack -- which is
        # the very TypeError this contract exists to prevent, raised from the wrong place.
        with mock.patch.object(perception_module, "get_lig_mol", return_value=(None, 0)):
            with self.assertRaises(ValueError) as ctx:
                get_tmc_mol(_FIXTURE, 0, with_stereo=False)
        self.assertIn("get_lig_mol failed", str(ctx.exception))

    @unittest.skipUnless(_FIXTURE.exists(), f"fixture missing: {_FIXTURE}")
    def test_broken_fixture_perceives_a_degenerate_graph_under_stable_metal_ac(self):
        """Documents the OIN_STABLE_METAL_AC consequence on this broken geometry.

        Not an endorsement -- a pin. If a future sanity gate makes this raise again, or the
        capping order changes, this test fails and the module docstring gets revisited rather
        than the behaviour drifting unobserved.
        """
        from rdkit import Chem

        result = get_tmc_mol(_FIXTURE, 0, with_stereo=False)
        mol = result[0] if isinstance(result, tuple) else result
        frags = Chem.GetMolFrags(mol)
        self.assertGreater(
            len(frags),
            1,
            "expected the degenerate multi-fragment perception this fixture now produces",
        )
        bare_protons = sum(
            1 for f in frags if len(f) == 1 and mol.GetAtomWithIdx(f[0]).GetAtomicNum() == 1
        )
        self.assertGreater(
            bare_protons,
            0,
            "the documented symptom is stranded bare-proton fragments; if that is gone, "
            "perception improved and this pin (and the module docstring) should be updated",
        )


if __name__ == "__main__":
    unittest.main()
