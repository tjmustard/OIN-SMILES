"""build_contract_mol must honor the OIN as the source of truth for sp3 stereo.

The round-trip re-encodes the generated structure through the *contract mol*
(``get_oin_string(gen_result.mol, coords)`` -- what the dataset harness does), NOT
``XYZToSMILES.convert``. ``build_contract_mol`` perceives chirality from the embed
geometry via ``AssignStereochemistryFrom3D``, which stamps a tetrahedral tag on
*every* chiral-looking centre in that one arbitrary conformer -- including centres
the OIN left UNspecified. The fully-sanitised forward encoder (``get_tmc_mol``)
drops those, so without a guard the re-encode invents an ``@`` the input never
emitted and the round trip fails on a centre that was never specified (KAPCEM:
0 in -> 4 out).

These tests pin that an OIN-unspecified sp3 centre stays unspecified through the
contract-mol re-encode, while a specified one is preserved.
"""

import unittest

import numpy as np
from rdkit import Chem

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.utils.xyz2mol import get_oin_string

# KAPCEM_comp_0's own forward encode: a Co complex whose tetradentate
# bis(hydroxylamine)diamine ligand carries NO sp3 '@' (the centres are left
# unspecified). The generated single conformer freezes those carbons into a
# chiral-looking geometry, so AssignStereochemistryFrom3D tags all four.
OIN_UNSPECIFIED = "[Co_SPY].c1ccc(P{0}(c2ccccc2)c2ccccc2)cc1.CC(N{1}O)C(C)N{3}CCCN{2}C(C)C(C)N{4}O"

# A square-planar Pt with a genuinely-specified free stereocentre three bonds from
# the donor (the existing sp3-carry fixture). Must round-trip WITH its '@'.
OIN_SPECIFIED = "[Pt_SPL].C[C@H](O)CN{0}.[Cl]{1}.[Cl]{2}.N{3}"

# AHEBEV_comp_0: a benzylic carbon on an eta6-arene bound to Cr. Its CIP label
# FLIPS between the metal-present contract mol (R) and the metal-free fragment (S),
# so the generator's metal-present flip loop mis-orients it. recover() must
# re-orient it on the metal-free fragment against the template rdCIPLabeler label,
# so the two specified centres round-trip with matching @ (canonical key equal).
OIN_ETA_ADJACENT = (
    "[Cr_TET].Cc{0}1[cH]{0}[cH]{0<}[cH]{0}[cH]{0}c{0}1[C@H]([C@H](O)C(C)(C)C)"
    "N(O)Cc1ccccc1.C{1}#O.C{2}#O.C{3}#O"
)

# BABWAD_comp_0: a menthyl carbon bonded to an eta5-cyclopentadienyl. rdCIPLabeler
# gives OPPOSITE R/S for that carbon depending on whether the Cp is left aromatic
# vs charged/kekulized, so the orientation stamp must be taken on the aromatic
# form (the form the fragment recover() re-orients against emits). Guards
# _template_sp3_label's aromatic-preserving sanitize.
OIN_CP_ADJACENT = (
    "[Ti_TET].CC(C)[C@@H]1CC[C@@H](C)C[C@H]1c{0}1[cH]{0}[cH]{0>}[cH]{0}[cH]{0}1."
    "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1c{1}1[cH]{1}[cH]{1>}[cH]{1}[cH]{1}1.[Cl]{2}.[Cl]{3}"
)


def _reencode_via_contract(oin):
    """Generate 3D then re-encode through the contract mol, as the harness does."""
    result = OIN3DGenerator(engine="metallogen", optimizer="ff").generate(oin)
    lines = result.xyz.splitlines()
    n = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + n)])
    assert result.mol is not None, "expected a non-None contract mol for this OIN"
    return get_oin_string(Chem.Mol(result.mol), coords)


class TestContractMolHonorsOINSpecification(unittest.TestCase):
    def test_unspecified_sp3_is_not_invented(self):
        """An OIN with no sp3 '@' must not gain one through the contract re-encode."""
        oin2 = _reencode_via_contract(OIN_UNSPECIFIED)
        self.assertNotIn("@", oin2, f"contract re-encode invented sp3 stereo: {oin2}")

    def test_unspecified_sp3_is_deterministic(self):
        """The seeded embed must give the same (still-unspecified) re-encode twice."""
        first = _reencode_via_contract(OIN_UNSPECIFIED)
        second = _reencode_via_contract(OIN_UNSPECIFIED)
        self.assertEqual(first, second)
        self.assertNotIn("@", first)

    def test_specified_sp3_is_preserved(self):
        """A genuinely specified free stereocentre must keep its '@'."""
        oin2 = _reencode_via_contract(OIN_SPECIFIED)
        self.assertIn("@", oin2, f"contract re-encode dropped specified stereo: {oin2}")

    def test_eta_adjacent_sp3_round_trips(self):
        """A metal/eta-adjacent specified sp3 centre must round-trip (correct @)."""
        oin2 = _reencode_via_contract(OIN_ETA_ADJACENT)
        self.assertEqual(
            canonical_roundtrip_key(OIN_ETA_ADJACENT),
            canonical_roundtrip_key(oin2),
            f"eta-adjacent sp3 stereo mis-oriented: {oin2}",
        )

    def test_cp_adjacent_sp3_round_trips(self):
        """A carbon bonded to an aromatic eta-Cp must round-trip (aromatic-form CIP)."""
        oin2 = _reencode_via_contract(OIN_CP_ADJACENT)
        self.assertEqual(
            canonical_roundtrip_key(OIN_CP_ADJACENT),
            canonical_roundtrip_key(oin2),
            f"Cp-adjacent sp3 stereo mis-oriented: {oin2}",
        )


if __name__ == "__main__":
    unittest.main()
