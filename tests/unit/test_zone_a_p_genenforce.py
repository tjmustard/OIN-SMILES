"""Stereo Phase 4 (Zone-A P Stereocenter Encoding) — generation-side tests.

Covers MiniPRD_ZoneA_P_GenEnforce.md ("MiniPRD-B") Sec.5 Tests 1-6 (US-B1,
US-B2, US-B3). Test 7 (Q3 molassembler investigation) is a Candidate-Artifact
finding recorded in spec/worklog/NOTES.md, not a unit test.

Fixture: tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz -- RhCl2(DIPAMP), both P atoms
genuine CIP stereocentres, single bidentate fragment (both P atoms coreside
in the SAME fragment -- this is what makes the "co-resident stereocenter
retains its configuration" assertion in Test 3 meaningful).
"""

import os
import sys
import time
import unittest
import warnings
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles import OINStereoWarning, XYZToSMILES
from oinsmiles.core.chirality import _metal_present_cip_label
from oinsmiles.generation import molassembler_adapter as ma
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.oin_parser import OINParser, OINVector, ParsedOIN

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_DIPAMP_XYZ = os.path.join(_FIXTURES_DIR, "Rh-RR-DIPAMP-Cl2.xyz")

# The DIPAMP ligand fragment, chirality markers stripped, used to locate the
# fragment-local atom index of each P atom (same numbering `_stitch_fragment`
# sees, since `[P@]`/`[P@@]`/`P` all parse to exactly one SMILES atom).
_DIPAMP_FRAG_NO_TAGS = "Cc1ccccc1P(CCP(c1ccccc1)c1ccccc1C)c1ccccc1"


def _dipamp_frag_p_local_idxs() -> list[int]:
    mol = Chem.MolFromSmiles(_DIPAMP_FRAG_NO_TAGS, sanitize=False)
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]


def _p_cip_codes_by_idx(mol: Chem.Mol) -> dict:
    """Metal-present rdCIPLabeler CIP codes for every P atom, keyed by index.

    Uses the shared `_metal_present_cip_label` diagnostic helper
    (core/chirality.py) rather than a raw `AssignStereochemistryFrom3D` call:
    the metal<->P bond is DATIVE in BOTH `get_tmc_mol()` output and the
    adapter's assembled `combined_mol`, but the donor/acceptor DIRECTION
    differs between the two (P is donor in `get_tmc_mol()`; the adapter's
    `_assemble_combined_mol` adds the bond metal->P, making the metal the
    donor) -- raw from-3D perception is silently sensitive to that direction
    (a pre-existing, documented gap; see `_metal_present_cip_label`'s
    docstring), so this test helper goes through the SAME apples-to-apples
    recipe used everywhere else in this MiniPRD for both sides of any
    comparison.
    """
    out = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            out[atom.GetIdx()] = _metal_present_cip_label(mol, atom.GetIdx())
    return out


class TestZoneAPParsedOINPassthrough(unittest.TestCase):
    """Task 1 (audit, not edit): [P@]{0} / [P@@]{1>} survive intact at the
    ParsedOIN level -- MiniPRD-A's Test 8 already covers the lower
    (parse_inline_string) level; this asserts the SAME property one layer up,
    at OINParser.parse()'s output, which is what `_template_generate`
    actually consumes.
    """

    def test_bracket_p_chiral_tags_survive_to_parsed_oin(self):
        inline = "[Rh_SPL].c1ccccc1[P@]{0}(C)c1ccccc1.[P@@]{1}.[Cl]{2}.[Cl]{3}"
        parsed = OINParser().parse(inline)

        self.assertIn("[P@]", parsed.fragments[1])
        self.assertEqual(parsed.fragments[2], "[P@@]")

        # atom_in_fragment_idx for the tagged P in fragment 1 must point at
        # the P atom itself (index 6 in "c1ccccc1[P@](C)c1ccccc1": 6 ring
        # atoms at indices 0-5, P at index 6).
        p_vector = next(v for v in parsed.vectors if v.fragment_idx == 1)
        self.assertEqual(p_vector.atom_in_fragment_idx, 6)
        frag_mol = Chem.MolFromSmiles(parsed.fragments[1], sanitize=False)
        self.assertEqual(frag_mol.GetAtomWithIdx(6).GetAtomicNum(), 15)
        self.assertNotEqual(
            frag_mol.GetAtomWithIdx(6).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED
        )


class TestZoneAPEnantiomerDiscrimination(unittest.TestCase):
    """Test 1 (US-B1.1/B1.2): DIPAMP OIN (fixed seed) regenerates a
    self-consistent metal-present CIP pair on both P atoms; flipping @<->@@
    on both P atoms in the OIN inverts both regenerated labels.
    """

    def test_flip_both_p_tags_inverts_both_regenerated_cip_labels(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)

        self.assertRegex(oin1, r"\[P@@?\]\{0\}")
        self.assertRegex(oin1, r"\[P@@?\]\{1\}")

        oin1_flipped = (
            oin1.replace("[P@@]", "\0TMP\0").replace("[P@]", "[P@@]").replace("\0TMP\0", "[P@]")
        )
        self.assertNotEqual(oin1, oin1_flipped)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            struct_a = OIN3DGenerator().generate(oin1)
            struct_b = OIN3DGenerator().generate(oin1_flipped)

        self.assertIsNotNone(struct_a.mol, "expected an assembled mol for struct_a")
        self.assertIsNotNone(struct_b.mol, "expected an assembled mol for struct_b")

        cip_a = sorted(_p_cip_codes_by_idx(struct_a.mol).values())
        cip_b = sorted(_p_cip_codes_by_idx(struct_b.mol).values())

        self.assertTrue(all(cip_a), f"struct_a missing CIP labels: {cip_a}")
        self.assertTrue(all(cip_b), f"struct_b missing CIP labels: {cip_b}")
        self.assertNotEqual(
            cip_a, cip_b, f"flipping both P tags did not invert regenerated CIP: {cip_a} vs {cip_b}"
        )


class TestZoneAPLosslessRoundTrip(unittest.TestCase):
    """Test 2 (US-B2): regenerated-structure CIP-from-3D matches the
    original fixture's CIP-from-3D.

    NOTE on scope: MiniPRD-B Test 2 also names full OIN-string byte-stability
    as an acceptance criterion. That is NOT achievable today: regenerating
    DIPAMP triggers PRE-EXISTING, unrelated generation-fidelity bugs already
    documented in test_stereo_roundtrip_diagnostics.py::test_p_stereocenter_roundtrip
    (geo_code drift SPL->SPY, a spurious C-C->C=C bond-order artifact, and a
    missing-Cl/extra-H perception issue on the re-derived topology) -- none of
    which are Zone-A-P-specific, and fixing them is out of MiniPRD-B's scope
    (it touches xyz2mol's bond-order perception on regenerated coordinates,
    not the adapter's stereo enforcement). That test remains
    `@unittest.expectedFailure`, documenting the gap honestly. This test
    covers the criterion that IS in scope and IS satisfied: the CIP oracle
    on the regenerated metal-present complex matches the original.
    """

    def test_regenerated_metal_present_cip_matches_original(self):
        from oinsmiles.utils.xyz2mol import get_tmc_mol

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            original_mol, _xyz = get_tmc_mol(_DIPAMP_XYZ, 0, with_stereo=False)
            Chem.SanitizeMol(original_mol)
            original_cip = sorted(_p_cip_codes_by_idx(original_mol).values())

            structure = OIN3DGenerator().generate(oin1)

        self.assertIsNotNone(structure.mol, "expected an assembled mol")
        regenerated_cip = sorted(_p_cip_codes_by_idx(structure.mol).values())

        self.assertTrue(all(original_cip), f"original fixture missing CIP labels: {original_cip}")
        self.assertTrue(all(regenerated_cip), f"regenerated missing CIP labels: {regenerated_cip}")
        self.assertEqual(
            original_cip,
            regenerated_cip,
            f"regenerated metal-present CIP does not match the original: "
            f"{original_cip} vs {regenerated_cip}",
        )


class TestZoneAPForcedMisEmbedCorrection(unittest.TestCase):
    """Test 3 (US-B1.3): via the Task-7 injection point (`_stitch_fragment`'s
    `_test_flip_chiral_idx`), force a mis-embed of exactly ONE Zone-A P atom
    on the FIRST placement attempt only. Enforcement must correct it within
    <=3 re-embed attempts, and the co-resident P atom (SAME bidentate
    fragment) must retain its configuration throughout -- proof no
    mirror/improper transform was applied (SuperPRD B2/B3).
    """

    def test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)
            baseline_structure = OIN3DGenerator().generate(oin1)
        self.assertIsNotNone(baseline_structure.mol)
        baseline_cip = _p_cip_codes_by_idx(baseline_structure.mol)

        p_local_idxs = _dipamp_frag_p_local_idxs()
        self.assertEqual(len(p_local_idxs), 2)
        target_p_local_idx = p_local_idxs[0]

        orig_stitch = ma._stitch_fragment
        call_count = {"n": 0}

        def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
            call_count["n"] += 1
            # Force the mis-embed ONLY on the very first call to the P-P
            # bidentate fragment (the injection point never fires again,
            # simulating a one-off ETKDG error rather than an unsatisfiable tag).
            if call_count["n"] == 1 and target_p_local_idx in binding_idxs:
                kwargs["_test_flip_chiral_idx"] = target_p_local_idx
            return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

        with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator().generate(oin1)

        self.assertIsNotNone(structure.mol, "expected an assembled mol after enforcement")
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            stereo_warnings,
            [],
            f"enforcement should have SUCCEEDED (no persistent-mismatch warning): "
            f"{[str(w.message) for w in stereo_warnings]}",
        )

        final_cip = _p_cip_codes_by_idx(structure.mol)
        self.assertEqual(
            final_cip,
            baseline_cip,
            "enforcement must converge back to the SAME CIP labels as the "
            f"unforced baseline (no residual mis-embed, no over-correction): "
            f"baseline={baseline_cip} final={final_cip}",
        )
        # At least one retry attempt must have actually happened (i.e. the
        # forced mismatch was real, not a no-op).
        self.assertGreater(call_count["n"], 3, "expected at least one re-embed attempt")


class TestZoneAPBoundedFailure(unittest.TestCase):
    """Test 4 (US-B1.4): an injection that flips the SAME atom on EVERY
    call (never satisfiable within the fixed-seed retry budget) must
    complete -- structure emitted, exactly one OINStereoWarning naming the
    offending atom, wall-clock far below the 60s ProcessPoolExecutor budget.
    """

    def test_persistent_mismatch_warns_once_and_completes_quickly(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)

        p_local_idxs = _dipamp_frag_p_local_idxs()
        target_p_local_idx = p_local_idxs[0]

        orig_stitch = ma._stitch_fragment

        def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
            if target_p_local_idx in binding_idxs:
                kwargs["_test_flip_chiral_idx"] = target_p_local_idx
            return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

        start = time.time()
        with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator().generate(oin1)
        elapsed = time.time() - start

        self.assertIsNotNone(
            structure.mol, "structure must still be emitted on persistent mismatch"
        )
        self.assertLess(elapsed, 30.0, "must complete far below the 60s ProcessPoolExecutor budget")

        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            len(stereo_warnings),
            1,
            f"expected exactly one OINStereoWarning, got: {[str(w.message) for w in stereo_warnings]}",
        )
        self.assertIn("could not be enforced", str(stereo_warnings[0].message))


class TestZoneAPFallbackObservability(unittest.TestCase):
    """Test 5 (US-B1.5, RISK-9): paths with no assembled RDKit mol skip
    enforcement + warn -- never a silent gap.
    """

    def test_template_path_with_no_assembled_mol_warns(self):
        """Simulates the eta-fallback case (`_template_generate` returns a
        result tuple whose mol is None) on an OIN that carries a Zone-A P
        tag -- `MolassemblerAdapter.generate()` must warn.
        """
        # Fragment must be a GENUINE CIP stereocentre (three chemically
        # distinct substituents: phenyl / ethyl / methyl) -- a symmetric P
        # (e.g. two identical phenyls) gets no expected label at all
        # (`_zone_a_p_expected_labels` correctly returns [] for it, the same
        # "silent, not a bug" outcome as BDPP/BDNN's negative controls) and
        # would never reach the fallback-warning code path under test here.
        parsed = ParsedOIN(
            smiles="[Rh].c1ccccc1[P@](CC)C.[Cl].[Cl]",
            fragments=["[Rh]", "c1ccccc1[P@](CC)C", "[Cl]", "[Cl]"],
            metal_fragment_idx=0,
            vectors=[
                OINVector(atom_idx=-1, vector=(1, 0, 0), fragment_idx=1, atom_in_fragment_idx=6),
                OINVector(atom_idx=-1, vector=(-1, 0, 0), fragment_idx=2, atom_in_fragment_idx=0),
                OINVector(atom_idx=-1, vector=(0, 1, 0), fragment_idx=3, atom_in_fragment_idx=0),
            ],
            original_oin="[Rh_TPL].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}",
            geo_code="TPL",
        )

        fake_result = ("3\ncomment\nRh 0 0 0\n", None, [])
        with mock.patch.object(ma, "_template_generate", return_value=fake_result):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = ma.MolassemblerAdapter().generate(parsed)

        self.assertIsNone(structure.mol)
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(len(stereo_warnings), 1)
        self.assertIn("stereo unenforced on fallback path", str(stereo_warnings[0].message))
        self.assertIn("atom 6", str(stereo_warnings[0].message))  # local P atom idx

    def test_dg_fallback_path_warns_when_zone_a_p_tag_present(self):
        """DG-path (geo_code='NON') never runs the assembled-complex
        enforcement -- must warn whenever the OIN carries a Zone-A P tag,
        independent of whether `_reconstruct_mol_from_smiles_and_xyz`
        happens to succeed.
        """
        from concurrent.futures import TimeoutError as FuturesTimeout  # noqa: F401
        from unittest.mock import MagicMock

        parsed = ParsedOIN(
            smiles="[Rh].c1ccccc1[P@](CC)C",
            fragments=["[Rh]", "c1ccccc1[P@](CC)C"],
            metal_fragment_idx=0,
            vectors=[
                OINVector(atom_idx=-1, vector=(1, 0, 0), fragment_idx=1, atom_in_fragment_idx=6),
            ],
            original_oin="[Rh_NON].c1ccccc1[P@]{0}(CC)C",
            geo_code="NON",
        )

        with mock.patch(
            "oinsmiles.generation.molassembler_adapter.ProcessPoolExecutor"
        ) as mock_pool_cls:
            mock_executor = MagicMock()
            mock_pool_cls.return_value.__enter__ = MagicMock(return_value=mock_executor)
            mock_pool_cls.return_value.__exit__ = MagicMock(return_value=False)
            mock_future = MagicMock()
            mock_future.result.return_value = {"error": "n/a", "ok": False}
            mock_executor.submit.return_value = mock_future

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with self.assertRaises(RuntimeError):
                    ma.MolassemblerAdapter().generate(parsed)

        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(len(stereo_warnings), 1)
        self.assertIn("Molassembler DG fallback path", str(stereo_warnings[0].message))


class TestZoneAPNoRegression(unittest.TestCase):
    """Test 6 (no regression): baseline generation suite unaffected;
    achiral/carbon-chirality generation paths and clean fixtures stay clean
    under -W error::OINStereoWarning.
    """

    def _assert_generates_clean(self, oin_string):
        with warnings.catch_warnings():
            warnings.simplefilter("error", OINStereoWarning)
            structure = OIN3DGenerator().generate(oin_string)
        self.assertTrue(len(structure.xyz) > 0)

    def test_cisplatin_generation_is_clean(self):
        self._assert_generates_clean("[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

    def test_ferrocene_generation_is_clean(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "ferrocene.xyz"))
        self._assert_generates_clean(oin)

    def test_dipamp_generation_is_clean(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(_DIPAMP_XYZ)
        self._assert_generates_clean(oin)

    def test_bdpp_generation_unaffected(self):
        """BDPP: P atoms are NOT CIP stereocentres (no _OIN_CIPCode_LP) --
        chirality lives on backbone carbon, an ordinary (pre-existing,
        untouched) path through `_stitch_fragment`."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDPP.xyz"))
        self._assert_generates_clean(oin)


if __name__ == "__main__":
    unittest.main()
