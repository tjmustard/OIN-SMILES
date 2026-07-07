"""Stereo Phase 4 (Zone-A P Stereocenter Encoding) tests.

Covers MiniPRD_ZoneA_P_Encode.md Sec.5 Tests 1, 2, 4, 5, 6, 7, 8, 9. Test 3
(negative controls) lives alongside the existing BDPP/BDNN goldens in
test_chiral_p.py / test_chiral_n.py.

Fixture: tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz -- RhCl2(DIPAMP), both P atoms
genuine CIP stereocentres (see file header comment for provenance).
"""

import hashlib
import os
import sys
import unittest
import warnings
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles import OINStereoWarning, XYZToSMILES
from oinsmiles.core import chirality as chirality_module
from oinsmiles.core.chirality import ChiralityRecoveryUtility, CIPAssigner
from oinsmiles.oin.inline import OINInlineHandler, SlotAssignment
from oinsmiles.utils.xyz2mol import get_tmc_mol

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_CANDIDATES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../candidate_outputs"))

_DIPAMP_XYZ = os.path.join(_FIXTURES_DIR, "Rh-RR-DIPAMP-Cl2.xyz")
_BDPP_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDPP.xyz")
_BDNN_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDNN.xyz")
_BINAP_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-R-BINAP.xyz")
_CISPLATIN_XYZ = os.path.join(_FIXTURES_DIR, "CisPlatin.xyz")
_FERROCENE_XYZ = os.path.join(_FIXTURES_DIR, "Ferrocene.xyz")

_BDPP_EXPECTED_OIN = (
    "[Pd_SPL].C[C@@H](C[C@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}"
)


def _dipamp_tmc_mol():
    """Fresh, sanitized, 3D full mol (metal present) for the DIPAMP fixture."""
    tmc_mol, _xyz = get_tmc_mol(_DIPAMP_XYZ, 0, with_stereo=False)
    Chem.SanitizeMol(tmc_mol)
    return tmc_mol


class TestZoneAPPropertyAssignment(unittest.TestCase):
    """Test 1: assign_all() stores _OIN_CIPCode_LP on both DIPAMP P atoms,
    idempotently."""

    def test_both_p_atoms_get_lp_cip_and_are_idempotent(self):
        mol = _dipamp_tmc_mol()
        p_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
        self.assertEqual(len(p_indices), 2, "DIPAMP fixture must have exactly 2 P atoms")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            CIPAssigner().assign_all(mol, diagnostics=False)

        first_pass = {}
        for idx in p_indices:
            label = mol.GetAtomWithIdx(idx).GetPropsAsDict().get("_OIN_CIPCode_LP")
            self.assertIn(label, ("R", "S"), f"atom {idx} missing/invalid _OIN_CIPCode_LP")
            first_pass[idx] = label

        # Idempotence: calling assign_all() again must yield identical properties.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            CIPAssigner().assign_all(mol, diagnostics=False)

        for idx in p_indices:
            second_label = mol.GetAtomWithIdx(idx).GetPropsAsDict().get("_OIN_CIPCode_LP")
            self.assertEqual(
                second_label,
                first_pass[idx],
                f"atom {idx}: _OIN_CIPCode_LP changed across repeated assign_all() calls",
            )


class TestZoneAPEmittedString(unittest.TestCase):
    """Test 2: XYZToSMILES().convert() on DIPAMP emits @/@@ on both P{0}/P{1}."""

    def test_dipamp_oin_has_tags_on_both_p_slots(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(_DIPAMP_XYZ)

        self.assertRegex(oin, r"\[P@@?\]\{0\}", f"P{{0}} untagged in OIN: {oin}")
        self.assertRegex(oin, r"\[P@@?\]\{1\}", f"P{{1}} untagged in OIN: {oin}")


class TestZoneAPRawParity(unittest.TestCase):
    """Test 4 (RISK-1, Q1): compare the fragment's raw trivalent rdCIPLabeler
    CIP (before recover()'s verify-and-flip) against the stored
    _OIN_CIPCode_LP. Documents divergence rather than silently hiding it."""

    def test_raw_fragment_cip_vs_stored_lp_label(self):
        """Capture the REAL fragment mol as the actual pipeline hands it to
        ``recover()`` (before any verify-and-flip runs) and compare its raw
        trivalent ``rdCIPLabeler`` CIP against the stored ``_OIN_CIPCode_LP``.

        Uses a monkeypatch on ``ChiralityRecoveryUtility.recover`` purely as
        a capture point -- the real ``recover()`` implementation still runs
        immediately afterwards, so the full pipeline behaviour (and its
        final OIN string) is completely unaffected.
        """
        from rdkit.Chem import rdCIPLabeler

        captured_fragments = []
        orig_recover = ChiralityRecoveryUtility.recover

        def _capturing_recover(self, mol):
            captured_fragments.append(Chem.Mol(mol) if mol is not None else None)
            return orig_recover(self, mol)

        with mock.patch.object(ChiralityRecoveryUtility, "recover", _capturing_recover):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                XYZToSMILES().convert(_DIPAMP_XYZ)

        divergences = []
        checked_atoms = set()
        for frag_mol in captured_fragments:
            if frag_mol is None:
                continue
            for atom in frag_mol.GetAtoms():
                if atom.GetAtomicNum() != 15 or not atom.HasProp("_OIN_CIPCode_LP"):
                    continue
                raw_frag = Chem.Mol(frag_mol)
                Chem.AssignStereochemistry(raw_frag, cleanIt=True, force=True)
                try:
                    rdCIPLabeler.AssignCIPLabels(raw_frag)
                except Exception:  # noqa: BLE001 - diagnostic only
                    continue
                raw_cip = raw_frag.GetAtomWithIdx(atom.GetIdx()).GetPropsAsDict().get("_CIPCode")
                stored = atom.GetPropsAsDict().get("_OIN_CIPCode_LP")
                checked_atoms.add(atom.GetIdx())
                if raw_cip != stored:
                    divergences.append((atom.GetIdx(), raw_cip, stored))

        self.assertEqual(len(checked_atoms), 2, "expected to check both DIPAMP P atoms")
        self.assertEqual(
            divergences,
            [],
            "Raw fragment trivalent CIP diverges from the stored lone-pair "
            "label (SuperPRD RISK-1/Q1). This is not necessarily a bug -- "
            "recover()'s verify-and-flip is designed to absorb exactly this "
            "-- but it must be investigated and the final OIN string "
            "(TestZoneAPEmittedString) re-checked before assuming this test "
            f"failure is benign: {divergences}",
        )


class TestZoneAPDegradation(unittest.TestCase):
    """Test 5 (RISK-7, B4): guarded degradation when the dummy-metal copy
    cannot be built (e.g. an eta-ligand aromatic ring stranded by bond
    removal in a CpM(PR3)-type complex).

    Engineering a REAL RDKit SanitizeMol failure for this shape is fragile
    and RDKit-version-dependent; the CONTRACT under test -- convert()
    completes, no property stored, output unchanged, exactly one
    OINStereoWarning naming the atom -- is deterministically exercised here
    by monkeypatching ``_build_dummy_metal_copy`` to simulate the failure on
    exactly one eligible Zone-A P (the BDPP fixture's P{0}), using the real
    BDPP fixture so "convert() completes" is a real, full pipeline run.
    """

    def test_dummy_copy_failure_degrades_gracefully(self):
        real_eligible = chirality_module._eligible_zone_a_p
        forced_idx = {}

        def _one_atom_only(mol):
            found = real_eligible(mol)
            restricted = found[:1]
            if restricted:
                forced_idx["idx"] = restricted[0]
            return restricted

        def _simulate_dummy_copy_failure(mol, p_idx):
            # Reproduces the real _build_dummy_metal_copy contract on
            # failure (guarded try/except -> warn + return None) without
            # depending on engineering a real, RDKit-version-fragile
            # SanitizeMol exception for this shape.
            warnings.warn(
                OINStereoWarning(
                    f"atom {p_idx}: dummy-metal copy construction failed "
                    "(simulated for test) -- Zone-A lone-pair CIP not "
                    "computed; degrades to today's clearing behaviour."
                ),
                stacklevel=2,
            )
            return None

        with (
            mock.patch.object(chirality_module, "_eligible_zone_a_p", side_effect=_one_atom_only),
            mock.patch.object(
                chirality_module,
                "_build_dummy_metal_copy",
                side_effect=_simulate_dummy_copy_failure,
            ),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            oin = XYZToSMILES().convert(_BDPP_XYZ)

        self.assertEqual(
            oin,
            _BDPP_EXPECTED_OIN,
            "convert() output must be byte-identical to pre-feature output "
            "when the dummy-metal copy cannot be built",
        )

        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            len(stereo_warnings),
            1,
            f"expected exactly one OINStereoWarning, got: {[str(w.message) for w in stereo_warnings]}",
        )
        self.assertIn(str(forced_idx["idx"]), str(stereo_warnings[0].message))


class TestZoneAPBridgingGuard(unittest.TestCase):
    """Test 6 (B7): a P bonded to >=2 metal atoms (bridging phosphide) gets
    no _OIN_CIPCode_LP and triggers exactly one OINStereoWarning; assign_all()
    completes without raising."""

    @staticmethod
    def _build_bridging_p_mol():
        rw = Chem.RWMol()
        p_idx = rw.AddAtom(Chem.Atom(15))
        rh1_idx = rw.AddAtom(Chem.Atom(45))
        rh2_idx = rw.AddAtom(Chem.Atom(45))
        c_idx = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(p_idx, rh1_idx, Chem.BondType.DATIVE)
        rw.AddBond(p_idx, rh2_idx, Chem.BondType.DATIVE)
        rw.AddBond(p_idx, c_idx, Chem.BondType.SINGLE)

        mol = rw.GetMol()
        conf = Chem.Conformer(mol.GetNumAtoms())
        coords = {
            p_idx: (0.0, 0.0, 0.0),
            rh1_idx: (1.8, 1.8, 1.8),
            rh2_idx: (1.8, -1.8, -1.8),
            c_idx: (-1.5, 1.5, -1.5),
        }
        for idx, xyz in coords.items():
            conf.SetAtomPosition(idx, xyz)
        mol.AddConformer(conf)
        Chem.SanitizeMol(mol)
        return mol, p_idx

    def test_bridging_p_gets_no_property_and_one_warning(self):
        mol, p_idx = self._build_bridging_p_mol()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = CIPAssigner().assign_all(mol, diagnostics=False)

        self.assertIsNotNone(result, "assign_all() must complete, not raise")
        self.assertFalse(
            mol.GetAtomWithIdx(p_idx).HasProp("_OIN_CIPCode_LP"),
            "bridging P must not receive a lone-pair CIP property",
        )
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(len(stereo_warnings), 1)
        self.assertIn(str(p_idx), str(stereo_warnings[0].message))


class TestZoneAPWarningGate(unittest.TestCase):
    """Test 7: clean fixtures pass under -W error::OINStereoWarning -- i.e.
    a real DeprecationWarning-style promotion of OINStereoWarning to an
    exception must NOT fire for any of these known-clean structures."""

    def _assert_clean(self, xyz_path):
        with warnings.catch_warnings():
            warnings.simplefilter("error", OINStereoWarning)
            oin = XYZToSMILES().convert(xyz_path)
        self.assertIsInstance(oin, str)
        self.assertGreater(len(oin), 0)

    def test_dipamp_is_clean(self):
        self._assert_clean(_DIPAMP_XYZ)

    def test_bdpp_is_clean(self):
        self._assert_clean(_BDPP_XYZ)

    def test_bdnn_is_clean(self):
        self._assert_clean(_BDNN_XYZ)

    def test_binap_is_clean(self):
        self._assert_clean(_BINAP_XYZ)

    def test_cisplatin_is_clean(self):
        self._assert_clean(_CISPLATIN_XYZ)

    def test_ferrocene_is_clean(self):
        self._assert_clean(_FERROCENE_XYZ)


class TestZoneAPParseAdjacency(unittest.TestCase):
    """Test 8 (B10): parse_inline_string on a string containing [P@]{0} and
    [P@@]{1>} -- slot, winding, and SMILES atom-index adjacency must all be
    correct. _count_smiles_atoms_before must treat a bracket-chirality token
    ([P@], [P@@]) as a single atom, same as any other bracket atom."""

    def test_bracket_p_chiral_tags_with_slot_and_winding(self):
        inline = "[Pt_SPL].[P@]{0}.[P@@]{1>}"
        smiles, geometry, vectors = OINInlineHandler.parse_inline_string(inline)

        self.assertEqual(geometry, "SPL")
        self.assertEqual(smiles, "[Pt].[P@].[P@@]")
        self.assertEqual(len(vectors), 2)
        self.assertEqual(vectors[0], SlotAssignment(1, 0, 0, None))
        self.assertEqual(vectors[1], SlotAssignment(2, 0, 1, ">"))


class TestZoneAPCandidateArtifact(unittest.TestCase):
    """Test 9 (Novel -- Candidate Artifact): the DIPAMP OIN string is written
    to tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt with a provenance
    line (fixture path + sha256) and is NOT silently promoted to
    tests/fixtures/ -- HITL review per SuperPRD Sec.9 is a separate,
    human-gated step this test cannot and must not perform."""

    def test_candidate_artifact_exists_with_provenance(self):
        candidate_path = os.path.join(_CANDIDATES_DIR, "Rh-RR-DIPAMP-Cl2_oin.txt")
        self.assertTrue(
            os.path.isfile(candidate_path),
            f"Candidate Artifact missing: {candidate_path}",
        )

        with open(candidate_path) as f:
            content = f.read()

        with open(_DIPAMP_XYZ, "rb") as f:
            expected_sha256 = hashlib.sha256(f.read()).hexdigest()

        self.assertIn(
            "tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz",
            content,
            "provenance line must name the canonical fixture path",
        )
        self.assertIn(
            expected_sha256,
            content,
            "provenance line must include the fixture's current sha256",
        )
        self.assertIn(
            "PENDING",
            content,
            "HITL sign-off status must be explicitly marked pending until a "
            "human records a verdict in spec/worklog/",
        )

        first_line = content.splitlines()[0]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            current_oin = XYZToSMILES().convert(_DIPAMP_XYZ)
        self.assertEqual(
            first_line,
            current_oin,
            "candidate file's OIN string must match the current pipeline output "
            "(stale candidate artifact -- regenerate via Task 13's script)",
        )

    def test_candidate_artifact_not_a_trusted_fixture(self):
        """The duplicate integration fixture is deleted (Task 12); the
        candidate file must live only under candidate_outputs/, never under
        fixtures/, until HITL sign-off promotes it."""
        self.assertFalse(
            os.path.isfile(
                os.path.join(os.path.dirname(__file__), "../integration/Rh-RR-DIPAMP-Cl2.xyz")
            ),
            "duplicate integration fixture should have been deleted (Task 12)",
        )


if __name__ == "__main__":
    unittest.main()
