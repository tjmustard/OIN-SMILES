"""Stereo Phase 4 (Zone-A P Stereocenter Encoding) — generation-side tests.

Covers MiniPRD_ZoneA_P_GenEnforce.md ("MiniPRD-B") Sec.5 Tests 1-6 (US-B1,
US-B2, US-B3). Test 7 (Q3 molassembler investigation) is a Candidate-Artifact
finding recorded in spec/worklog/NOTES.md, not a unit test.

Fixture: tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz -- RhCl2(DIPAMP), both P atoms
genuine CIP stereocentres, single bidentate fragment (both P atoms coreside
in the SAME fragment). Still used by Tests 1/2 (enantiomer discrimination,
lossless round-trip), which only exercise CIP-oracle agreement and do not
depend on which placement path (template/Kabsch vs Molassembler DG) the
adapter chooses.

TASK-31 update (2026-07-03): DIPAMP's bite is incompatible with the
Kabsch/template placement path, so the adapter now correctly routes it to
the Molassembler DG fallback (see `_stitch_fragment`'s non-binding-H guard
in molassembler_adapter.py). The DG path never runs the assembled-complex
Zone-A P enforcement machinery (`_verify_zone_a_p`), so DIPAMP no longer
exercises Tests 3/4/6 below. TASK-32 retargets those to SYNTHETIC,
hand-written inline OINs built around a MONODENTATE P-stereocenter
(`c1ccccc1[P@](...)C` on a `_TET` geometry) that verifiably stays on the
template path -- see `_MONO_P_*` fixtures below. Test 3's "co-resident
stereocenter not disturbed" check now uses a single ligand carrying BOTH a
Zone-A P stereocenter and a directly-bonded carbon stereocenter
(`[P@]([C@@H](C)CC)C`), in place of DIPAMP's second P atom.
"""

import os
import sys
import time
import unittest
import warnings
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from oinsmiles import OINStereoWarning, XYZToSMILES
from oinsmiles.core.chirality import _lp_cip_label, _metal_present_cip_label
from oinsmiles.generation import molassembler_adapter as ma
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.oin_parser import OINParser, OINVector, ParsedOIN

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_INTEGRATION_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../integration"))
_DIPAMP_XYZ = os.path.join(_FIXTURES_DIR, "Rh-RR-DIPAMP-Cl2.xyz")

# MiniPRD-C (Zone-A P SPL Dummy-Metal Embed) fixtures.
# Monodentate 3D round-trip fixture (Test 4, ACCEPTANCE/HARD GATE): built
# independently in Avogadro (not derived from any oinsmiles output), a
# genuine P stereocentre (o-tolyl / phenyl / methyl / metal) on Rh_SPL with
# a PMe3 co-ligand. Lives in tests/integration/ only (not dual-copied to
# tests/fixtures/ -- an existing, pre-MiniPRD-C inconsistency in this repo's
# fixture layout, not introduced here).
_MONO_P_SPL_REAL_XYZ = os.path.join(_FIXTURES_DIR, "Rh-Single-Chiral-Phosphine.xyz")

# Two incompatible-bite bidentate 3D fixtures (Test 8): both confirmed
# empirically (this MiniPRD's own spike) to route to the DG fallback via
# ee0b3f0's bite-distortion guard -- the MiniPRD's Task 8 text calls the
# first "(compatible-bite)", which is incorrect (see spec/process/ for the
# resolution): DIPAMP is the incompatible-bite case ee0b3f0 was written for.
# Per user decision, both fixtures are used as incompatible-bite regression
# cases; the compatible-bite bidentate 3D case remains a known gap.
_DIPAMP_PH_XYZ = os.path.join(_FIXTURES_DIR, "Rh-RR-DIPAMP-Ph-Cl2.xyz")

# Same monodentate P-stereocentre fragment as _MONO_P_FRAG_NO_TAGS below, but
# on [Pt_SPL] instead of [Ni_TET] -- exercises the CONFIRMED SPL bug this
# MiniPRD fixes directly (self-oracle, no XYZ needed).
_MONO_P_SPL_OIN = "[Pt_SPL].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}.[Cl]{3}"
_MONO_P_SPL_FRAG_TAGGED_R = "c1ccccc1[P@](CC)C"
_MONO_P_SPL_FRAG_TAGGED_S = "c1ccccc1[P@@](CC)C"
_MONO_P_CORESIDENT_SPL_OIN = "[Pt_SPL].c1ccccc1[P@]{0}([C@@H](C)CC)C.[Cl]{1}.[Cl]{2}.[Cl]{3}"

# The DIPAMP ligand fragment, chirality markers stripped, used to locate the
# fragment-local atom index of each P atom (same numbering `_stitch_fragment`
# sees, since `[P@]`/`[P@@]`/`P` all parse to exactly one SMILES atom).
_DIPAMP_FRAG_NO_TAGS = "Cc1ccccc1P(CCP(c1ccccc1)c1ccccc1C)c1ccccc1"


def _dipamp_frag_p_local_idxs() -> list[int]:
    mol = Chem.MolFromSmiles(_DIPAMP_FRAG_NO_TAGS, sanitize=False)
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]


# --- TASK-32: synthetic monodentate-P fixtures that stay on the template
# (Kabsch-enforcing) path -- verified empirically (2026-07-03) to build a
# mol with 0 OINStereoWarnings at baseline. Denticity-independent per the
# TASK-32 worklog: `_verify_zone_a_p` runs at the assembled-complex level in
# `_template_generate` for ANY fragment carrying a Zone-A P tag that reaches
# the template path, so a monodentate P exercises the same machinery DIPAMP
# used to.

# Simple monodentate P-stereocenter: phenyl / ethyl / methyl on P (plus the
# metal as the 4th, dative substituent). Used by Tests 4/6 (persistent
# mismatch, no-regression).
_MONO_P_OIN = "[Ni_TET].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}.[Cl]{3}"
_MONO_P_FRAG_NO_TAGS = "c1ccccc1P(CC)C"

# Same P-stereocenter, but one substituent is itself a carbon stereocenter
# (-CH(CH3)CH2CH3) directly bonded to P -- gives a single monodentate ligand
# with BOTH a Zone-A P stereocenter AND a co-resident (non-metal-binding)
# carbon stereocenter, replacing DIPAMP's second P atom for Test 3.
_MONO_P_CORESIDENT_OIN = "[Ni_TET].c1ccccc1[P@]{0}([C@@H](C)CC)C.[Cl]{1}.[Cl]{2}.[Cl]{3}"
_MONO_P_CORESIDENT_FRAG_NO_TAGS = "c1ccccc1P(C(C)CC)C"


def _mono_p_frag_p_local_idx() -> int:
    mol = Chem.MolFromSmiles(_MONO_P_FRAG_NO_TAGS, sanitize=False)
    idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    assert len(idxs) == 1
    return idxs[0]


def _mono_p_coresident_frag_p_local_idx() -> int:
    mol = Chem.MolFromSmiles(_MONO_P_CORESIDENT_FRAG_NO_TAGS, sanitize=False)
    idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    assert len(idxs) == 1
    return idxs[0]


def _p_and_co_resident_c_cip_codes(mol: Chem.Mol) -> dict:
    """CIP codes for a Zone-A P stereocenter AND its directly-bonded
    co-resident carbon stereocenter, keyed by ``("P"|"C", atom_idx)``.

    The P side uses the same metal-present diagnostic recipe as
    `_p_cip_codes_by_idx` (the metal<->P dative bond direction differs
    between `get_tmc_mol()` output and the adapter's assembled
    `combined_mol`, so raw from-3D perception on P is direction-sensitive --
    see `_metal_present_cip_label`'s docstring). The co-resident carbon is
    NOT bonded to the metal, so it has no such ambiguity: a direct
    `rdCIPLabeler` recompute on the assembled mol is unambiguous.
    """
    out: dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15 and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            out[("P", atom.GetIdx())] = _metal_present_cip_label(mol, atom.GetIdx())
            for nbr in atom.GetNeighbors():
                if (
                    nbr.GetAtomicNum() == 6
                    and nbr.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                ):
                    mol_copy = Chem.Mol(mol)
                    rdCIPLabeler.AssignCIPLabels(mol_copy)
                    out[("C", nbr.GetIdx())] = (
                        mol_copy.GetAtomWithIdx(nbr.GetIdx()).GetPropsAsDict().get("_CIPCode")
                    )
    return out


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
            struct_a = OIN3DGenerator(engine="legacy").generate(oin1)
            struct_b = OIN3DGenerator(engine="legacy").generate(oin1_flipped)

        self.assertIsNotNone(struct_a.mol, "expected an assembled mol for struct_a")
        self.assertIsNotNone(struct_b.mol, "expected an assembled mol for struct_b")

        cip_a = sorted(_p_cip_codes_by_idx(struct_a.mol).values())
        cip_b = sorted(_p_cip_codes_by_idx(struct_b.mol).values())

        self.assertTrue(all(cip_a), f"struct_a missing CIP labels: {cip_a}")
        self.assertTrue(all(cip_b), f"struct_b missing CIP labels: {cip_b}")
        self.assertNotEqual(
            cip_a, cip_b, f"flipping both P tags did not invert regenerated CIP: {cip_a} vs {cip_b}"
        )


def _single_p_global_idx(mol: Chem.Mol) -> int:
    """Return the GetIdx() of the sole P atom in *mol* (fixtures below have exactly one)."""
    p_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    assert len(p_idxs) == 1, f"expected exactly one P atom, found {len(p_idxs)}"
    return p_idxs[0]


class TestZoneAPSPLDummyEmbed(unittest.TestCase):
    """MiniPRD-C Test 1 (ACCEPTANCE): the CONFIRMED SPL bug this MiniPRD
    fixes -- spec/worklog/SPL-P-enforcement-decision.md. Before the
    dummy-metal embed, SPL's metal-present CIP was fixed by which face the
    metal happened to land on, not by the ETKDG re-embed seed, so only ONE
    of [P@]/[P@@] ever enforced cleanly (the other always warned "could not
    be enforced"). Both tags must now enforce cleanly, with like-for-like
    (lone-pair) labels matching the OIN-encoded intent and OPPOSITE
    metal-present CIPs (enantiomer discrimination).
    """

    def test_both_spl_tags_enforce_cleanly_with_opposite_metal_present_cip(self):
        metal_present_labels = {}
        for tag, frag in (
            ("[P@]", _MONO_P_SPL_FRAG_TAGGED_R),
            ("[P@@]", _MONO_P_SPL_FRAG_TAGGED_S),
        ):
            oin = f"[Pt_SPL].c1ccccc1{tag}{{0}}(CC)C.[Cl]{{1}}.[Cl]{{2}}.[Cl]{{3}}"
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator(engine="legacy").generate(oin)
            self.assertIsNotNone(structure.mol, f"expected an assembled mol for {tag}")

            stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
            enforce_warnings = [
                w for w in stereo_warnings if "could not be enforced" in str(w.message)
            ]
            self.assertEqual(
                enforce_warnings,
                [],
                f"{tag} must enforce cleanly on SPL (this MiniPRD's own bug fix): "
                f"{[str(w.message) for w in enforce_warnings]}",
            )

            p_global_idx = _single_p_global_idx(structure.mol)
            expected_lp = dict(ma._zone_a_p_expected_labels(frag))[_mono_p_frag_p_local_idx()]
            dummy = ma._build_dummy_metal_copy(structure.mol, p_global_idx)
            self.assertIsNotNone(dummy, f"dummy-metal copy failed for {tag}")
            measured_lp = _lp_cip_label(dummy, p_global_idx)
            self.assertEqual(
                measured_lp,
                expected_lp,
                f"{tag}: like-for-like LP label mismatch (measured={measured_lp}, "
                f"expected={expected_lp})",
            )

            metal_present_labels[tag] = _metal_present_cip_label(structure.mol, p_global_idx)

        self.assertTrue(all(metal_present_labels.values()), metal_present_labels)
        self.assertNotEqual(
            metal_present_labels["[P@]"],
            metal_present_labels["[P@@]"],
            f"opposite tags must give opposite metal-present CIPs: {metal_present_labels}",
        )

    def test_first_embed_enforces_without_entering_reembed_loop(self):
        """Test 6: with Tasks 1-6 correct, `_verify_zone_a_p` returns no
        mismatches on the FIRST embed for both tags on SPL, so the re-embed
        loop's attempt counter stays 0 -- correctness moved to embed time,
        the loop is a pure safety net (Task 7).
        """
        for tag in ("[P@]", "[P@@]"):
            oin = f"[Pt_SPL].c1ccccc1{tag}{{0}}(CC)C.[Cl]{{1}}.[Cl]{{2}}.[Cl]{{3}}"
            call_count = {"n": 0}
            orig_stitch = ma._stitch_fragment

            def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
                call_count["n"] += 1
                return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

            with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    structure = OIN3DGenerator(engine="legacy").generate(oin)
            self.assertIsNotNone(structure.mol)
            # Exactly one _stitch_fragment call for the single ligand
            # fragment (Cl atoms are single-atom fragments too, but the
            # P-bearing ligand is what matters): no retry attempt fired.
            self.assertEqual(
                call_count["n"],
                4,
                f"expected no re-embed retries for {tag} (4 fragments, 1 call each): "
                f"got {call_count['n']} _stitch_fragment calls",
            )


class TestZoneAPSPLForcedMisEmbedCorrection(unittest.TestCase):
    """MiniPRD-C Test 3 (co-resident safety on SPL) and Test 7
    (loop-with-dummy): extends MiniPRD-B's forced-mis-embed test to SPL,
    and explicitly asserts no Z=0/`*` atom survives the re-embed loop's own
    dummy attach+strip (Task 6's postcondition, exercised on the loop path).
    """

    def test_single_atom_mis_embed_on_spl_is_corrected_without_mirroring_co_resident(self):
        oin1 = _MONO_P_CORESIDENT_SPL_OIN
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            baseline_structure = OIN3DGenerator(engine="legacy").generate(oin1)
        self.assertIsNotNone(baseline_structure.mol)
        baseline_cip = _p_and_co_resident_c_cip_codes(baseline_structure.mol)
        self.assertTrue(any(k[0] == "P" for k in baseline_cip))
        self.assertTrue(any(k[0] == "C" for k in baseline_cip))

        target_p_local_idx = _mono_p_coresident_frag_p_local_idx()
        orig_stitch = ma._stitch_fragment
        call_count = {"n": 0}

        def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
            call_count["n"] += 1
            if call_count["n"] == 1 and target_p_local_idx in binding_idxs:
                kwargs["_test_flip_chiral_idx"] = target_p_local_idx
            return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

        with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator(engine="legacy").generate(oin1)

        self.assertIsNotNone(structure.mol, "expected an assembled mol after enforcement")
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            stereo_warnings,
            [],
            f"enforcement should have SUCCEEDED: {[str(w.message) for w in stereo_warnings]}",
        )

        final_cip = _p_and_co_resident_c_cip_codes(structure.mol)
        self.assertEqual(
            final_cip,
            baseline_cip,
            f"enforcement must converge to the SAME CIP labels: baseline={baseline_cip} "
            f"final={final_cip}",
        )
        self.assertGreater(call_count["n"], 4, "expected at least one re-embed attempt")

        # Task 6 postcondition, exercised on the LOOP path (Test 7): no Z=0
        # dummy leaked into the final assembled mol, and atom count is
        # exactly what it would be with no dummy ever having existed.
        self.assertFalse(
            any(a.GetAtomicNum() == 0 for a in structure.mol.GetAtoms()),
            "no Z=0 dummy atom may survive the re-embed loop",
        )
        n_declared = int(structure.xyz.strip().split("\n")[0])
        self.assertNotIn(
            "*",
            {line.split()[0] for line in structure.xyz.strip().split("\n")[2 : 2 + n_declared]},
            "no '*' (Z=0) token may appear in the written XYZ block",
        )


class TestZoneAPByteStableRoundTrip(unittest.TestCase):
    """MiniPRD-C Test 4 (ACCEPTANCE, HARD GATE, C5): XYZ->OIN->XYZ->OIN on
    the monodentate real-3D fixture -> second OIN byte-identical to the
    first. Per C5, the fixture's absolute configuration is trusted as-is
    (RISK-C3 accepted); reviewer sign-off recorded post-hoc.
    """

    def test_rh_single_chiral_phosphine_round_trips_byte_identically(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            oin1 = XYZToSMILES().convert(_MONO_P_SPL_REAL_XYZ)

        with warnings.catch_warnings(record=True) as caught2:
            warnings.simplefilter("always")
            structure = OIN3DGenerator(engine="legacy").generate(oin1)

        self.assertIsNotNone(structure.mol, "expected an assembled mol")
        stereo_warnings = [
            w for w in list(caught) + list(caught2) if issubclass(w.category, OINStereoWarning)
        ]
        self.assertEqual(
            stereo_warnings,
            [],
            f"expected clean generation: {[str(w.message) for w in stereo_warnings]}",
        )

        import tempfile

        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            f.write(structure.xyz)
            tmp_path = f.name
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                oin2 = XYZToSMILES().convert(tmp_path)
        finally:
            os.unlink(tmp_path)

        self.assertEqual(
            oin1, oin2, f"round trip not byte-identical:\n  oin1={oin1}\n  oin2={oin2}"
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

            structure = OIN3DGenerator(engine="legacy").generate(oin1)

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
    <=3 re-embed attempts, and a CO-RESIDENT stereocenter in the SAME
    fragment must retain its configuration throughout -- proof no
    mirror/improper transform was applied (SuperPRD B2/B3).

    TASK-32 retarget (2026-07-03): DIPAMP no longer reaches the template
    (Kabsch-enforcing) path (TASK-31 routes it to DG), so the co-resident
    stereocenter it used to supply (its second P atom) is replaced by
    `_MONO_P_CORESIDENT_OIN` -- a single monodentate phosphine ligand with
    BOTH a Zone-A P stereocenter AND a directly-bonded carbon stereocenter
    (`[P@]([C@@H](C)CC)C`), verified empirically to build a mol with 0
    OINStereoWarnings at baseline on the template path.
    """

    def test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident(self):
        oin1 = _MONO_P_CORESIDENT_OIN
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            baseline_structure = OIN3DGenerator(engine="legacy").generate(oin1)
        self.assertIsNotNone(baseline_structure.mol)
        baseline_cip = _p_and_co_resident_c_cip_codes(baseline_structure.mol)
        self.assertTrue(
            any(k[0] == "P" for k in baseline_cip), "expected a Zone-A P CIP label in baseline"
        )
        self.assertTrue(
            any(k[0] == "C" for k in baseline_cip),
            "expected a co-resident carbon CIP label in baseline",
        )

        target_p_local_idx = _mono_p_coresident_frag_p_local_idx()

        orig_stitch = ma._stitch_fragment
        call_count = {"n": 0}

        def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
            call_count["n"] += 1
            # Force the mis-embed ONLY on the very first call to the P
            # fragment (the injection point never fires again, simulating a
            # one-off ETKDG error rather than an unsatisfiable tag).
            if call_count["n"] == 1 and target_p_local_idx in binding_idxs:
                kwargs["_test_flip_chiral_idx"] = target_p_local_idx
            return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

        with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator(engine="legacy").generate(oin1)

        self.assertIsNotNone(structure.mol, "expected an assembled mol after enforcement")
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            stereo_warnings,
            [],
            f"enforcement should have SUCCEEDED (no persistent-mismatch warning): "
            f"{[str(w.message) for w in stereo_warnings]}",
        )

        final_cip = _p_and_co_resident_c_cip_codes(structure.mol)
        self.assertEqual(
            final_cip,
            baseline_cip,
            "enforcement must converge back to the SAME CIP labels (Zone-A P AND "
            "co-resident carbon) as the unforced baseline (no residual mis-embed, "
            f"no over-correction, no mirroring): baseline={baseline_cip} final={final_cip}",
        )
        # At least one retry attempt must have actually happened (i.e. the
        # forced mismatch was real, not a no-op).
        self.assertGreater(call_count["n"], 3, "expected at least one re-embed attempt")


class TestZoneAPBoundedFailure(unittest.TestCase):
    """Test 4 (US-B1.4): an injection that flips the SAME atom on EVERY
    call (never satisfiable within the fixed-seed retry budget) must
    complete -- structure emitted, exactly one OINStereoWarning naming the
    offending atom, wall-clock far below the 60s ProcessPoolExecutor budget.

    TASK-32 retarget (2026-07-03): vehicle changed from DIPAMP (now routes
    to DG, TASK-31) to `_MONO_P_OIN`, a synthetic monodentate P-stereocenter
    OIN verified to stay on the template/Kabsch path and enforce cleanly at
    baseline (0 warnings) -- so a persistent forced mismatch here still
    exercises the SAME assembled-complex `_verify_zone_a_p` bounded-retry
    machinery DIPAMP used to.
    """

    def test_persistent_mismatch_warns_once_and_completes_quickly(self):
        oin1 = _MONO_P_OIN
        target_p_local_idx = _mono_p_frag_p_local_idx()

        orig_stitch = ma._stitch_fragment

        def _wrapped(frag_smiles, binding_idxs, target_positions, **kwargs):
            if target_p_local_idx in binding_idxs:
                kwargs["_test_flip_chiral_idx"] = target_p_local_idx
            return orig_stitch(frag_smiles, binding_idxs, target_positions, **kwargs)

        start = time.time()
        with mock.patch.object(ma, "_stitch_fragment", side_effect=_wrapped):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                structure = OIN3DGenerator(engine="legacy").generate(oin1)
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

    def test_dg_fallback_total_failure_raises_without_stereo_warning(self):
        """MiniPRD-C (Task 5, C2) supersedes the old unconditional up-front
        warning: enforcement now runs AFTER a `bonded_mol` exists (a
        set-based `_zone_a_p_measured_labels_dg` comparison, since the
        connected SMILES is canonicalised and per-atom index tracking isn't
        available here). When Molassembler itself fails outright (``ok:
        False``), generation raises `RuntimeError` before any mol exists to
        verify -- the exception itself is the loud, non-silent signal
        (RISK-9's intent), so no separate OINStereoWarning fires here.
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
        self.assertEqual(stereo_warnings, [])


class TestZoneAPIncompatibleBiteDGEnforcement(unittest.TestCase):
    """MiniPRD-C Test 8 (C2): bidentate Zone-A-P chelates whose bite is
    incompatible with template placement (ee0b3f0 routes these to the DG
    fallback) must STILL be enforced -- Task 5's set-based verify-and-
    reseed loop on the DG-produced mol. Both fixtures here are confirmed
    (this MiniPRD's own spike) to be the incompatible-bite case; no
    compatible-bite bidentate 3D fixture currently exists (known gap, user
    decision 2026-07-03: not blocking this MiniPRD).
    """

    def _assert_enforces_cleanly_via_dg(self, xyz_path):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(xyz_path)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            structure = OIN3DGenerator(engine="legacy").generate(oin)

        self.assertTrue(len(structure.xyz) > 0)
        stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
        self.assertEqual(
            stereo_warnings,
            [],
            f"expected clean DG-path enforcement (no fallback/could-not-enforce "
            f"warning): {[str(w.message) for w in stereo_warnings]}",
        )
        self.assertIsNotNone(structure.mol, "expected a reconstructed mol on the DG path")
        self.assertFalse(
            any(a.GetAtomicNum() == 0 for a in structure.mol.GetAtoms()),
            "no Z=0 dummy atom may survive onto the DG-path mol",
        )
        n_declared = int(structure.xyz.strip().split("\n")[0])
        self.assertNotIn(
            "*",
            {line.split()[0] for line in structure.xyz.strip().split("\n")[2 : 2 + n_declared]},
            "no '*' (Z=0) token may appear in the written XYZ block",
        )

    def test_dipamp_enforces_cleanly_via_dg_fallback(self):
        self._assert_enforces_cleanly_via_dg(_DIPAMP_XYZ)

    def test_dipamp_ph_enforces_cleanly_via_dg_fallback(self):
        self._assert_enforces_cleanly_via_dg(_DIPAMP_PH_XYZ)


class TestZoneAPNoRegression(unittest.TestCase):
    """Test 6 (no regression): baseline generation suite unaffected;
    achiral/carbon-chirality generation paths and clean fixtures stay clean
    under -W error::OINStereoWarning.
    """

    def _assert_generates_clean(self, oin_string):
        with warnings.catch_warnings():
            warnings.simplefilter("error", OINStereoWarning)
            structure = OIN3DGenerator(engine="legacy").generate(oin_string)
        self.assertTrue(len(structure.xyz) > 0)

    def test_cisplatin_generation_is_clean(self):
        self._assert_generates_clean("[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

    def test_ferrocene_generation_is_clean(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "Ferrocene.xyz"))
        self._assert_generates_clean(oin)

    def test_monodentate_p_generation_is_clean(self):
        """TASK-32 retarget (2026-07-03): renamed from
        `test_dipamp_generation_is_clean`. DIPAMP now correctly routes to
        the DG fallback (TASK-31) and no longer exercises the
        template/Kabsch-path enforcement machinery -- see
        `test_dipamp_dg_fallback_warns_honestly` below for DIPAMP's own
        (now honest, non-silent) behaviour. This test keeps the original
        "a P-stereocenter complex generates clean on the enforcing template
        path" intent using `_MONO_P_OIN`, a synthetic monodentate
        P-stereocenter fixture verified to stay on the template path.
        """
        self._assert_generates_clean(_MONO_P_OIN)

    def test_bdpp_generation_unaffected(self):
        """BDPP: P atoms are NOT CIP stereocentres (no _OIN_CIPCode_LP) --
        chirality lives on backbone carbon, an ordinary (pre-existing,
        untouched) path through `_stitch_fragment`."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            oin = XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDPP.xyz"))
        self._assert_generates_clean(oin)

    def test_tag_free_goldens_never_enter_dummy_embed_branch(self):
        """MiniPRD-C Task 10 (C4.4): tag-free goldens stay byte-identical
        AND the dummy-embed branch is proven NEVER ENTERED for them via a
        branch-entry spy on `_attach_dummy_metal` -- not merely an output
        byte-diff (a byte-diff alone couldn't distinguish "branch never
        ran" from "branch ran and happened to cancel out").
        """
        golden_oins = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            golden_oins.append(("cisplatin", "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"))
            golden_oins.append(
                (
                    "transplatin",
                    XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "TransPlatin.xyz")),
                )
            )
            golden_oins.append(
                (
                    "cis_ptcl2en",
                    XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "Cis-PtCl2(en).xyz")),
                )
            )
            golden_oins.append(
                ("ferrocene", XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "Ferrocene.xyz")))
            )
            golden_oins.append(
                (
                    "fac_irppy3",
                    XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "fac-Ir(ppy)3.xyz")),
                )
            )
            golden_oins.append(
                (
                    "mer_irppy3",
                    XYZToSMILES().convert(os.path.join(_FIXTURES_DIR, "mer-Ir(ppy)3.xyz")),
                )
            )

        from oinsmiles.core import chirality as chirality_mod

        call_count = {"n": 0}
        orig_attach = chirality_mod._attach_dummy_metal

        def _spy(mol, p_idx):
            call_count["n"] += 1
            return orig_attach(mol, p_idx)

        for name, oin in golden_oins:
            call_count["n"] = 0
            with mock.patch.object(ma, "_attach_dummy_metal", side_effect=_spy):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    structure = OIN3DGenerator(engine="legacy").generate(oin)
            stereo_warnings = [w for w in caught if issubclass(w.category, OINStereoWarning)]
            self.assertEqual(stereo_warnings, [], f"{name}: unexpected stereo warning(s)")
            self.assertTrue(len(structure.xyz) > 0, f"{name}: expected non-empty XYZ")
            self.assertEqual(
                call_count["n"],
                0,
                f"{name}: dummy-embed branch was entered ({call_count['n']} calls) for a "
                "tag-free golden -- must be gated strictly on "
                "_zone_a_p_expected_labels(frag_smiles) non-empty",
            )


if __name__ == "__main__":
    unittest.main()
