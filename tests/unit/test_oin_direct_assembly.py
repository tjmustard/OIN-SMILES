"""SL2 oin-direct-winding: direct OIN->MetalComplex assembly (Part A).

Covers the assembly plumbing that bypasses the winding-lossy m-SMILES bridge:
 - the ``convert_parsed_to_msmiles`` refactor stays byte-identical (m-SMILES path),
 - ``om.get_om_from_parsed`` builds a MetalComplex equivalent to the m-SMILES path
   for a non-haptic complex (metal-first ordering, same slots/charge/multiplicity),
 - the direct path generates a non-haptic complex key-equal to the default path,
 - the ``oin_direct`` flag is read from ff_params and OIN_DIRECT_ASSEMBLY.
"""

import os
import unittest
from unittest import mock

import numpy as np
from rdkit import Chem

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.metallogen_adapter import (
    DEFAULT_SELECT_POOL,
    MetalloGenAdapter,
    _prepare_ligand_fragments,
    convert_parsed_to_msmiles,
)
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import om
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.utils.xyz2mol import get_oin_string

CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
FERROCENE = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"


def _reencode(result):
    lines = result.xyz.splitlines()
    n = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + n)])
    assert result.mol is not None, "expected a non-None contract mol"
    return get_oin_string(Chem.Mol(result.mol), coords)


class TestMSmilesByteIdentity(unittest.TestCase):
    """The _prepare_ligand_fragments refactor must not change the m-SMILES output."""

    def test_ferrocene_msmiles_unchanged(self):
        self.assertEqual(
            convert_parsed_to_msmiles(OINParser().parse(FERROCENE)),
            "[Fe]|[cH:1]1[cH:1][cH-:1][cH:1][cH:1]1|[cH:2]1[cH:2][cH-:2][cH:2][cH:2]1|2_linear",
        )

    def test_cisplatin_msmiles_stays_cis(self):
        # cis = the two Cl on adjacent (not opposite) slots; a stable regression pin.
        ms = convert_parsed_to_msmiles(OINParser().parse(CISPLATIN))
        self.assertTrue(ms.startswith("[Pt]|"))
        self.assertTrue(ms.endswith("|4_square_planar"))


class TestDirectConstructorEquivalence(unittest.TestCase):
    """get_om_from_parsed matches get_om_from_modified_smiles for a non-haptic complex."""

    def test_cisplatin_complex_equivalent(self):
        parsed = OINParser().parse(CISPLATIN)
        metal_frag, specs, geo = _prepare_ligand_fragments(parsed)
        direct = om.get_om_from_parsed(metal_frag, specs, geo)
        ref = om.get_om_from_modified_smiles(convert_parsed_to_msmiles(parsed))

        self.assertEqual(direct.metal_index, 0, "metal-first invariant")
        self.assertEqual(direct.center_atom.get_element(), ref.center_atom.get_element())
        self.assertEqual(direct.chg, ref.chg)
        self.assertEqual(direct.multiplicity, ref.multiplicity)
        self.assertEqual(len(direct.ligands), len(ref.ligands))
        self.assertEqual(
            [bi[1] for lig in direct.ligands for bi in lig.binding_infos],
            [bi[1] for lig in ref.ligands for bi in lig.binding_infos],
            "same donor -> slot assignment",
        )
        # Non-haptic: no winding attached.
        self.assertTrue(all(lig.winding is None for lig in direct.ligands))

    def test_direct_generation_key_equal_to_default(self):
        default = OIN3DGenerator(engine="metallogen", optimizer="ff").generate(CISPLATIN)
        direct = OIN3DGenerator(
            engine="metallogen", optimizer="ff", ff_params={"oin_direct": True}
        ).generate(CISPLATIN)
        self.assertEqual(
            canonical_roundtrip_key(_reencode(default)),
            canonical_roundtrip_key(_reencode(direct)),
        )


class TestFlagGating(unittest.TestCase):
    """oin_direct is off by default; ff_params and OIN_DIRECT_ASSEMBLY both enable it."""

    def test_ff_params_flag(self):
        self.assertTrue(MetalloGenAdapter(ff_params={"oin_direct": True})._oin_direct_enabled())
        self.assertFalse(MetalloGenAdapter()._oin_direct_enabled())

    def test_env_flag(self):
        with mock.patch.dict(os.environ, {"OIN_DIRECT_ASSEMBLY": "1"}):
            self.assertTrue(MetalloGenAdapter()._oin_direct_enabled())
        with mock.patch.dict(os.environ, {"OIN_DIRECT_ASSEMBLY": "0"}):
            self.assertFalse(MetalloGenAdapter()._oin_direct_enabled())

    def test_direct_dg_on_by_default(self):
        # v0.4.4: OIN-direct assembly feeds the DG embed as the DEFAULT path (an A/B matched
        # the m-SMILES path byte-for-byte per molecule). Guard a silent revert to off.
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("OIN_DIRECT_DG", None)
            self.assertTrue(MetalloGenAdapter()._direct_dg_enabled())
        with mock.patch.dict(os.environ, {"OIN_DIRECT_DG": "0"}):
            self.assertFalse(MetalloGenAdapter()._direct_dg_enabled())
        self.assertFalse(MetalloGenAdapter(ff_params={"direct_dg": False})._direct_dg_enabled())


class TestHapticPoolCollapse(unittest.TestCase):
    """Under oin_direct, a haptic complex routes to rigid placement with NO widening."""

    def test_haptic_collapses_pool_and_routes_kabsch(self):
        import oinsmiles.generation.metallogen_adapter as mad

        parsed = OINParser().parse(FERROCENE)
        captured = {}

        class _Stop(Exception):
            pass

        def fake(*args, **kwargs):
            captured["args"] = args
            captured["kwargs"] = kwargs
            raise _Stop()

        adapter = MetalloGenAdapter(optimizer="ff", ff_params={"oin_direct": True})
        with mock.patch.object(mad, "generate_3d_structures", fake):
            with self.assertRaises(_Stop):
                adapter.generate(parsed)

        kw = captured["kwargs"]
        # Pool collapsed to the default (5), NOT the ETA_SELECT_POOL=16 search width;
        # UFF pre-pool not doubled to 2*pool.
        self.assertEqual(kw["num_conformers"], DEFAULT_SELECT_POOL)
        self.assertLessEqual(kw["uff_pool_size"], 10)
        # Direct assembly (pre-built complex, no m-SMILES string) and haptic->kabsch.
        self.assertIsNone(captured["args"][0])
        self.assertIsNotNone(kw["metal_complex"])
        self.assertTrue(kw["ff_params"].get("kabsch_only"))

    def test_default_path_widens_pool_for_contrast(self):
        # The off-flag path DOES widen for the winding search -- proof the collapse
        # above is the oin_direct behaviour, not a no-op.
        import oinsmiles.generation.metallogen_adapter as mad

        parsed = OINParser().parse(FERROCENE)
        captured = {}

        class _Stop(Exception):
            pass

        def fake(*args, **kwargs):
            captured["kwargs"] = kwargs
            raise _Stop()

        adapter = MetalloGenAdapter(optimizer="ff")
        with mock.patch.object(mad, "generate_3d_structures", fake):
            with self.assertRaises(_Stop):
                adapter.generate(parsed)
        self.assertGreater(captured["kwargs"]["num_conformers"], DEFAULT_SELECT_POOL)


if __name__ == "__main__":
    unittest.main()
