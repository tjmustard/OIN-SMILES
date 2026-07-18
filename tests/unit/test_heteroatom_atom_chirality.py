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

# KAGXUM_comp_0: a methine carbon bonded to a substituted Cp AND a FUSED fluorenyl.
# For a fused haptic ring the _oin_fragment_templates object's RemoveHs corrupts the
# aromatic state, flipping rdCIPLabeler even under the aromatic-preserving sanitize,
# so the label must be taken on a FRESH re-parse of the template SMILES. Guards
# _template_sp3_label's re-parse.
OIN_FUSED_RING_ADJACENT = (
    "[Zr_TET].Cc{0}1[cH]{0>}[cH]{0}c{0}(C)c{0}1[C@H](c1ccccc1)"
    "c{1>}1c{1}2ccccc{1}2c{1}2ccccc{1}12.[Cl]{2}.[Cl]{3}"
)

# GUXPIA_comp_0: a Zone-A phosphorus DONOR whose lone-pair CIP label is
# representation-sensitive (bonded to an aromatic phenyl + a Cp arm). The lone-pair
# stamp (_template_lp_label) and recover()'s Zone-A P branch both read it on the same
# aromatic-preserving re-parse, so the donor round-trips at the correct handedness.
OIN_P_DONOR_AROMATIC = (
    "[Ni_SPL].[CH]1[CH][CH]C([P@@]{0}(Cc2ccc(CP{2}(c3ccccc3)c3ccccc3)n{1}2)"
    "c2ccccc2)[CH][CH]1.S{3}c1ccccc1"
)

# HUGSEI_comp_0: a Mo carbonyl with two thioether S DONORS carrying no stereo. The
# generator's AssignStereochemistryFrom3D stamps a spurious high-coordination
# permutation tag ([S@SP3]/[S@SP1]/[S@TB9H]) on each donor S from the metal-present
# geometry -- a tag the input crystal geometry never produced. recover() must clear it:
# a thioether S is not a genuine stereocentre once the metal is removed. (Fix A.)
OIN_DONOR_S_THIOETHER = "[Mo_OCT].CCS{2}CC[NH]{4}CCS{0}CC.C{1}#O.C{3}#O.C{5}#O"

# REPZUJ_comp_0: a GENUINE sulfonimidoyl S(VI) stereocentre (S bonded to aryl, =O, a
# metal-donor O, and =N). The donor O becomes a radical [O] in the metal-free fragment,
# which rdCIPLabeler refuses to rank -- so recover()'s _SP3_CIP_PROP re-orientation
# returned no label and was silently skipped, leaving the arbitrary embed handedness.
# H-filling the open valence makes the CIP computable, in the same convention the
# template stamp reads, so the centre round-trips at the correct handedness. (Fix B.)
OIN_SULFONIMIDOYL_S = (
    "[Rh_TPY].Cc1ccc([S@](O{0})(=O)=N{2}C[C@H](c2ccccc2)c2ccccc2-c2ccccn{3}2)cc1."
    "Cc{1}1c{1>}(C)c{1}(C)c{1}(C)c{1}1C"
)

# POYJIX_comp_0: a GENUINE quaternary ammonium N+ (bridged bicyclic) is a real
# tetrahedral stereocentre. The forward encoder keeps it (CIPAssigner stores
# _OIN_CIPCode), but build_contract_mol did not stamp N, so recover()'s 4-neighbour
# no-CIP fallback cleared it ([N@@+] -> [N+]). build_contract_mol now routes a degree-4
# N+ through the same metal-free _SP3_CIP_PROP re-orientation as C/Si/S. (Fix B.)
OIN_QUATERNARY_N = (
    "[Ni_SPL].COc1ccccc1N{0}C(=O)CN1C{1}N(C)C=C1.COc1ccccc1N{3}C(=O)C[N@@+]12C{2}=N(C)C1[CH-]2"
)

# XILZID_comp_0: an sp3 lactone carbon bonded to a metal-bound alkene DONOR (a
# valence-deficient carbon in the fragment). Its template CIP and fragment CIP diverged
# unless the donor's open valence is normalised identically on both sides, so the
# _SP3_CIP_PROP re-orientation mis-fired and inverted the centre. Both the stamp and the
# recover() comparison now read the label through the same fill-first reparse. (Fix B.)
OIN_ALKENE_DONOR_C = (
    "[Au_LIN].CC(C)c1cc(C(C)C)c(-c2ccccc2P{0}(C2CCCCC2)C2CCCCC2)c(C(C)C)c1."
    "CC1=C{1}[C@H](C(C)C)OC1=O"
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

    def test_fused_ring_adjacent_sp3_round_trips(self):
        """A carbon bonded to a fused eta-fluorenyl must round-trip (re-parsed CIP)."""
        oin2 = _reencode_via_contract(OIN_FUSED_RING_ADJACENT)
        self.assertEqual(
            canonical_roundtrip_key(OIN_FUSED_RING_ADJACENT),
            canonical_roundtrip_key(oin2),
            f"fused-ring-adjacent sp3 stereo mis-oriented: {oin2}",
        )

    def test_aromatic_p_donor_lone_pair_round_trips(self):
        """A Zone-A P donor with an aromatic arm must round-trip (re-parsed LP CIP)."""
        oin2 = _reencode_via_contract(OIN_P_DONOR_AROMATIC)
        self.assertEqual(
            canonical_roundtrip_key(OIN_P_DONOR_AROMATIC),
            canonical_roundtrip_key(oin2),
            f"P-donor lone-pair mis-oriented: {oin2}",
        )

    def test_donor_thioether_s_stereo_is_not_invented(self):
        """A metal-donor thioether S must not gain a spurious high-coordination tag."""
        oin2 = _reencode_via_contract(OIN_DONOR_S_THIOETHER)
        self.assertNotIn("[S@", oin2, f"contract re-encode invented donor-S stereo: {oin2}")

    def test_sulfonimidoyl_s_round_trips(self):
        """A genuine sulfonimidoyl S(VI) stereocentre must round-trip (correct @)."""
        oin2 = _reencode_via_contract(OIN_SULFONIMIDOYL_S)
        self.assertEqual(
            canonical_roundtrip_key(OIN_SULFONIMIDOYL_S),
            canonical_roundtrip_key(oin2),
            f"sulfonimidoyl S mis-oriented: {oin2}",
        )

    def test_quaternary_ammonium_n_round_trips(self):
        """A genuine quaternary ammonium N+ stereocentre must round-trip (not cleared)."""
        oin2 = _reencode_via_contract(OIN_QUATERNARY_N)
        self.assertEqual(
            canonical_roundtrip_key(OIN_QUATERNARY_N),
            canonical_roundtrip_key(oin2),
            f"quaternary N+ dropped or mis-oriented: {oin2}",
        )

    def test_alkene_donor_carbon_round_trips(self):
        """An sp3 C bonded to a metal-bound alkene donor must round-trip (fill-first CIP)."""
        oin2 = _reencode_via_contract(OIN_ALKENE_DONOR_C)
        self.assertEqual(
            canonical_roundtrip_key(OIN_ALKENE_DONOR_C),
            canonical_roundtrip_key(oin2),
            f"alkene-donor-adjacent sp3 carbon mis-oriented: {oin2}",
        )


if __name__ == "__main__":
    unittest.main()
