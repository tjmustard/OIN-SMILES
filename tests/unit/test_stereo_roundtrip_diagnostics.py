"""Diagnostic round-trip tests: expose currently-silent stereo loss in the
OIN→XYZ (generation) direction.

The XYZ→OIN direction encodes ligand P/N stereocenters as ``@``/``@@`` in the
fragment SMILES and eta-ring winding as ``{n>}``/``{n<}`` slot markers. The
OIN→XYZ direction (``OIN3DGenerator`` in ``src/oinsmiles/generation/engine.py``)
currently does not round-trip these: winding direction is discarded by
``SLOT_REGEX`` (``src/oinsmiles/oin/inline.py:44`` only captures the numeric
slot rank, not the ``>``/``<`` suffix), P/N CIP codes are never read during
generation, and haptic ring FACE is chosen only by flipping the ring normal
toward the metal (``src/oinsmiles/generation/molassembler_adapter.py:627``).

These three tests measure exactly that gap: XYZ → OIN(1) → generate 3D →
XYZ(2) → OIN(2), and check whether OIN(1) == OIN(2). They are marked
``@unittest.expectedFailure`` so the suite stays green while documenting the
gap; they are the acceptance tests for ROADMAP-stereo.md Phases 1-4 (each
phase flips one of these to a hard assert once fixed).

This test module MEASURES ONLY — it must never modify anything under src/.
"""

import os
import re
import sys
import tempfile
import unittest

import numpy as np
import rdkit
import scipy
from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_CANDIDATES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../candidate_outputs"))

_BDPP_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDPP.xyz")
_BDNN_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDNN.xyz")
_FERROCENE_OIN_PATH = os.path.join(_CANDIDATES_DIR, "ferrocene_oin.txt")
_DIPAMP_XYZ = os.path.join(_FIXTURES_DIR, "Rh-RR-DIPAMP-Cl2.xyz")

# Stereo Phase 3 (Haptic Face Correction) fixtures.
_HALIDE_FACE_XYZ = os.path.join(
    os.path.dirname(__file__), "../integration/Ferrocene-halide-face.xyz"
)
_HALIDE_FACE_GOLDEN_PATH = os.path.join(_CANDIDATES_DIR, "Ferrocene-halide-face_oin.txt")
_CHIRALITY_WITNESS_OIN_PATH = os.path.join(_CANDIDATES_DIR, "ChiralityWitnessRing_oin.txt")
_ANSA_BASELINE_XYZ_PATH = os.path.join(
    _CANDIDATES_DIR, "AnsaMetallocene_TiCat1_before_baseline.xyz"
)
# Natural (unmodified) OIN encoding of tests/integration/TiCat1.xyz, hardcoded
# here so this test suite does not depend on the XYZ->OIN encoder's exact
# output to exercise the generation-side US-004 conflict/no-regression path.
_ANSA_OIN_NATURAL = (
    "[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0}[cH]{0}[cH]{0<}1)"
    "c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}"
)
# Same string with ring0's marker flipped -- forces a genuine US-004 conflict
# (ring0 wants a flip; ring1 carries no marker at all, so it "does not").
_ANSA_OIN_CONFLICT = _ANSA_OIN_NATURAL.replace("{0<}", "{0>}", 1)
# Both rings given a marker that mismatches their natural embedding -- forces
# a genuine US-004 *coherent* case (both rings agree they want a flip).
_ANSA_OIN_COHERENT_FIRE = _ANSA_OIN_NATURAL.replace("{0<}", "{0>}", 1).replace(
    "[cH]{1}1", "[cH]{1<}1", 1
)

# Stereo Phase 3 golden was blessed against these exact tool versions (the
# entire diagnostic rests on ETKDG seed-42 embedding determinism, which is
# not guaranteed stable across rdkit/numpy/scipy releases). pyproject.toml
# only pins lower bounds (">="), matching this repo's existing convention,
# so the exact blessed versions are recorded here instead, per the MiniPRD's
# "record in the test/worklog if pinning in pyproject.toml isn't
# appropriate" fallback.
_BLESSED_RDKIT_VERSION = "2025.09.3"
_BLESSED_NUMPY_VERSION = "2.2.6"
_BLESSED_SCIPY_VERSION = "1.15.3"

# Matches a bracketed P atom carrying an explicit chirality flag, e.g.
# "[P@]" or "[P@@]". SMILES always requires brackets to carry @/@@, so
# this is a reliable (if not exhaustive) sanity probe.
_P_CHIRAL_MARKER = "[P@"

# Matches an OIN slot marker, capturing the optional winding direction.
_SLOT_MARKER_RE = re.compile(r"\{\d+([><])?\}")


def _fragment_signature_and_windings(fragment: str) -> tuple:
    """Return (marker-stripped fragment string, [winding chars found]).

    The marker-stripped string is a ring's content-based identity that
    survives both a re-encoder heading-atom change (R2) and a fragment-order
    swap -- it depends only on ring substituent identity/sequence, never on
    which specific atom happens to carry the visible marker.
    """
    windings = [m.group(1) for m in _SLOT_MARKER_RE.finditer(fragment) if m.group(1)]
    stripped = _SLOT_MARKER_RE.sub("", fragment)
    return stripped, windings


def _ring_winding_by_signature(oin_string: str) -> dict:
    """Map each ligand ring's content signature to its single winding char.

    Skips the metal fragment (index 0) and any fragment with no winding
    marker at all (zero-marker rings, monodentate ligands, …).
    """
    parts = oin_string.split(".")
    result = {}
    for frag in parts[1:]:
        sig, windings = _fragment_signature_and_windings(frag)
        if not windings:
            continue
        assert len(windings) == 1, (
            f"expected exactly one winding marker per ring, found {windings} in {frag!r}"
        )
        result[sig] = windings[0]
    return result


def _p_cip_codes(mol) -> list:
    """Return the sorted list of RDKit _CIPCode values found on P atoms in *mol*."""
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    codes = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus
            cip = atom.GetPropsAsDict().get("_CIPCode")
            if cip:
                codes.append(cip)
    return sorted(codes)


def _generate_and_reencode(oin_string: str) -> str:
    """Run OIN(1) -> 3D structure -> XYZ file -> OIN(2), returning OIN(2).

    Raises whatever OIN3DGenerator.generate() or XYZToSMILES().convert()
    raise; callers are expected to catch and self.fail() with the message
    so expectedFailure records the crash as the diagnostic result.
    """
    generator = OIN3DGenerator()
    structure = generator.generate(oin_string)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as tmp_file:
        tmp_file.write(structure.xyz)
        tmp_path = tmp_file.name

    try:
        return XYZToSMILES().convert(tmp_path)
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


class TestStereoRoundTripDiagnostics(unittest.TestCase):
    """Diagnostics for OIN->XYZ->OIN stereo stability (ROADMAP-stereo.md)."""

    # NOTE (2026-07-03): This test was originally written with
    # @unittest.expectedFailure per the roadmap's Phase 1/2/4 assumption that
    # OIN->XYZ generation drops P-stereocenter @/@@ tags. It UNEXPECTEDLY
    # PASSED (OIN(1) == OIN(2), verbatim) — see spec/worklog/NOTES.md Log
    # entry for 2026-07-03. Root cause: in the BDPP fixture the chirality is
    # carried by backbone CARBON atoms (P atoms are not CIP stereocentres —
    # both carry two identical phenyl groups), and the ligand fragment's
    # @/@@ tags are passed straight through SMILES embedding rather than
    # re-derived from 3D geometry, so nothing is lost for this fixture.
    # The decorator is intentionally OFF: this is now a hard regression
    # assert. The corresponding roadmap phase can be downgraded from "fix
    # loss" to "add hard test" for this pathway; true P-stereocenter (P atom
    # itself is a CIP center) coverage still needs a dedicated fixture.
    def test_chiral_p_roundtrip(self):
        """Phase 1/2/4 roadmap item — downgraded to hard assert 2026-07-03
        after unexpected pass (see NOTES.md).

        P-stereocenter @/@@ tags (here: on backbone C atoms) must survive
        OIN -> 3D generation -> re-encoding.
        """
        oin1 = XYZToSMILES().convert(_BDPP_XYZ)

        # Sanity check: today's XYZ->OIN direction DOES encode ligand chiral
        # tags (this assertion should pass; failure here would indicate the
        # XYZ->OIN pipeline itself regressed, not the generation gap).
        self.assertIn("@", oin1, "Expected @/@@ ligand chiral tags in OIN(1)")

        try:
            oin2 = _generate_and_reencode(oin1)
        except Exception as exc:  # noqa: BLE001 - diagnostic capture
            self.fail(f"generation crashed: {exc!r}")
            return

        self.assertEqual(
            oin1,
            oin2,
            f"OIN round-trip mismatch (P-stereocenters lost):\n  OIN(1): {oin1}\n  OIN(2): {oin2}",
        )

    # NOTE (2026-07-03): Same downgrade as test_chiral_p_roundtrip above —
    # this test was originally @unittest.expectedFailure and UNEXPECTEDLY
    # PASSED (OIN(1) == OIN(2), verbatim). See spec/worklog/NOTES.md Log
    # entry for 2026-07-03. Root cause identical: BDNN's chirality lives on
    # backbone C atoms (N atoms are tertiary amines with two identical
    # phenyl groups, not CIP stereocentres), and the @/@@ tags pass through
    # SMILES embedding unchanged. Decorator intentionally OFF; now a hard
    # regression assert.
    def test_chiral_n_roundtrip(self):
        """Phase 1/2/4 roadmap item — downgraded to hard assert 2026-07-03
        after unexpected pass (see NOTES.md).

        N-stereocenter @/@@ tags (here: on backbone C atoms) must survive
        OIN -> 3D generation -> re-encoding.
        """
        oin1 = XYZToSMILES().convert(_BDNN_XYZ)

        # Sanity check: today's XYZ->OIN direction DOES encode ligand chiral
        # tags (this assertion should pass).
        self.assertIn("@", oin1, "Expected @/@@ ligand chiral tags in OIN(1)")

        try:
            oin2 = _generate_and_reencode(oin1)
        except Exception as exc:  # noqa: BLE001 - diagnostic capture
            self.fail(f"generation crashed: {exc!r}")
            return

        self.assertEqual(
            oin1,
            oin2,
            f"OIN round-trip mismatch (N-stereocenters lost):\n  OIN(1): {oin1}\n  OIN(2): {oin2}",
        )

    # TASK-20 Phase 2 diagnostic (2026-07-03): closes the gap TASK-10 left
    # open -- BDPP/BDNN carry chirality on backbone CARBON atoms, not on the
    # P/N donor atom itself. Fixture: tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz, an
    # (R,R)-DIPAMP RhCl2 complex where both P atoms are genuine CIP
    # stereocentres (phenyl / 2-methoxyphenyl / ethylene-bridge / metal
    # substituents), built independently of the OIN pipeline (Avogadro), not
    # derived from any oinsmiles output.
    #
    # UPDATE (Stereo Phase 4, 2026-07-03): the ENCODING-side gap described
    # below is fixed by MiniPRD_ZoneA_P_Encode.md -- CIPAssigner.assign_all()
    # now stores a fragment-local (lone-pair convention) CIP label
    # (``_OIN_CIPCode_LP``) on each Zone-A P bonded to exactly one metal, via
    # a dummy-metal copy computed while 3D is still present; recover() keeps
    # and verify-and-flips that tag instead of clearing it. The SANITY
    # assertion below (OIN(1) carries @/@@ on the P atom) now PASSES. The
    # round-trip assertion is still expected to fail: OIN->XYZ generation
    # enforcement of the P tag (re-embed on mismatch) is MiniPRD-B
    # (MiniPRD_ZoneA_P_GenEnforce.md, not yet implemented) -- today's
    # regeneration also perturbs unrelated things (geometry code SPL->SPY,
    # spurious H atoms, a C=C bond artifact), so OIN(2) still differs from
    # OIN(1) for reasons well beyond just the P tag. Decorator stays ON.
    @unittest.expectedFailure
    def test_p_stereocenter_roundtrip(self):
        """P-stereocenter (donor atom itself) round-trip: XYZ -> OIN(1) ->
        generate 3D -> XYZ(2) -> OIN(2).

        Diagnostic result (2026-07-03, pre-Phase-4): failed at the SANITY
        assertion, not the round-trip assertion. OIN(1) never carried @/@@
        on the P atom at all: ``ChiralityRecoveryUtility.recover()``
        unconditionally cleared the chiral tag on any P/N atom with
        ``total_degree < 4`` in the POST-FRAGMENTATION ligand mol ("Zone A"
        atoms). Because OIN's ligand-fragment SMILES excludes the metal by
        construction, any P/N atom that binds the metal directly always has
        exactly 3 neighbours in the fragment -- so it was ALWAYS Zone A and
        ALWAYS stripped, regardless of fixture quality.

        Diagnostic result (2026-07-03, post-Phase-4): the SANITY assertion
        now passes -- see the class-level NOTE above. The round-trip
        assertion still fails: this is now purely a GENERATION-side
        (OIN->XYZ) gap, tracked by MiniPRD-B.
        """
        oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)

        self.assertIn(
            _P_CHIRAL_MARKER,
            oin1,
            "OIN(1) has no @/@@ on the P atom -- this is an XYZ->OIN "
            f"ENCODING gap (Zone A stripping), not a generation gap. "
            f"OIN(1): {oin1}",
        )

        try:
            oin2 = _generate_and_reencode(oin1)
        except Exception as exc:  # noqa: BLE001 - diagnostic capture
            self.fail(f"generation crashed: {exc!r}")
            return

        self.assertEqual(
            oin1,
            oin2,
            f"OIN round-trip mismatch (P-stereocenter lost):\n  OIN(1): {oin1}\n  OIN(2): {oin2}",
        )

    # NOTE (Stereo Phase 4, 2026-07-03): this test was originally
    # @unittest.expectedFailure per the pre-Phase-4 assumption that OIN(1)
    # never carries @/@@ on the P atom (same encoding gap as
    # test_p_stereocenter_roundtrip above). It UNEXPECTEDLY PASSED in full
    # (encode-side fix + today's ETKDG already building the correct pyramid
    # from the trivalent [P@] tag, per SuperPRD_StereoPhase4 spike 2 -- no
    # MiniPRD-B re-embed enforcement needed for THIS assertion to hold).
    # Decorator intentionally OFF: hard regression assert. MiniPRD-B may
    # still be needed for cases where ETKDG mis-embeds and enforcement must
    # correct it; this test only requires the flip to produce a detectable
    # difference, which today's un-enforced path already delivers.
    def test_p_stereocenter_flip_inverts_cip(self):
        """Flipping the P `@`/`@@` in the input OIN must produce generated
        3D structures with OPPOSITE CIP at the P atom (RDKit oracle on
        ``GeneratedStructure.mol``, per ROADMAP-stereo.md H-1).

        Diagnostic result (2026-07-03, pre-Phase-4): blocked at the same
        precondition as ``test_p_stereocenter_roundtrip`` -- OIN(1) carried
        no @/@@ on the P atom to flip in the first place, so no "twin" input
        could even be constructed.

        Diagnostic result (2026-07-03, post-Phase-4): downgraded to a hard
        assert -- see the NOTE above.
        """
        oin1 = XYZToSMILES().convert(_DIPAMP_XYZ)

        self.assertIn(
            _P_CHIRAL_MARKER,
            oin1,
            "Cannot construct a flipped twin: OIN(1) carries no @/@@ on the "
            f"P atom (same encoding gap as test_p_stereocenter_roundtrip). OIN(1): {oin1}",
        )

        if "[P@@]" in oin1:
            oin1_flipped = oin1.replace("[P@@]", "[P@]", 1)
        else:
            oin1_flipped = oin1.replace("[P@]", "[P@@]", 1)
        self.assertNotEqual(oin1, oin1_flipped, "Flip must actually change the input string")

        try:
            struct_a = OIN3DGenerator().generate(oin1)
            struct_b = OIN3DGenerator().generate(oin1_flipped)
        except Exception as exc:  # noqa: BLE001 - diagnostic capture
            self.fail(f"generation crashed: {exc!r}")
            return

        if struct_a.mol is not None and struct_b.mol is not None:
            cip_a = _p_cip_codes(struct_a.mol)
            cip_b = _p_cip_codes(struct_b.mol)
            self.assertTrue(cip_a, "No P _CIPCode assignable on generated structure A")
            self.assertTrue(cip_b, "No P _CIPCode assignable on generated structure B")
            self.assertNotEqual(
                cip_a,
                cip_b,
                f"Flipped-input P CIP codes did not differ: A={cip_a} B={cip_b}",
            )
        else:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as fa:
                fa.write(struct_a.xyz)
                path_a = fa.name
            with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as fb:
                fb.write(struct_b.xyz)
                path_b = fb.name
            try:
                oin2a = XYZToSMILES().convert(path_a)
                oin2b = XYZToSMILES().convert(path_b)
            finally:
                for p in (path_a, path_b):
                    try:
                        os.unlink(p)
                    except OSError:
                        pass

            self.assertNotEqual(
                oin2a,
                oin2b,
                f"Flipped-input re-encoded OIN P tags did not differ:\n  A: {oin2a}\n  B: {oin2b}",
            )

    @unittest.skip(
        "Symmetry-impossibility (US-003): plain ferrocene's Cp rings are "
        "substituent-symmetric (every ring atom is a plain [cH], identical "
        "local signature) -- winding is not a geometric observable for a "
        "symmetric ring, so no correction (fired or otherwise) can ever make "
        "this fixture's re-encoded winding marker track the input. The "
        "generation-side correction correctly detects this and reports a "
        "runtime 'no-op' (identity, no wasted rotation) for both rings -- see "
        "test_haptic_face_symmetric_ring_is_inert_noop. This test is kept "
        "(not deleted) to document the honest, non-fakeable outcome; the "
        "substituted Ferrocene-halide-face fixture is the real acceptance "
        "test for Phase 3 (see test_haptic_face_* below)."
    )
    def test_haptic_face_winding(self):
        """Fixed by Phase 3: eta-ring winding markers ({n>} vs {n<}) must
        survive OIN -> 3D generation -> re-encoding, and flipping one ring's
        winding in the input OIN must produce a distinguishable output.

        Diagnostic result (2026-07-03): winding markers do not survive
        generation at all (SLOT_REGEX only captures the numeric slot rank,
        discarding the direction suffix) — see NOTES.md Log entry.

        Diagnostic result (Stereo Phase 3): the underlying mechanism now
        exists and works (see the substituted-ring tests below) but ferrocene
        itself is symmetric, so this specific fixture can never distinguish
        the two windings -- skipped rather than faked.
        """
        with open(_FERROCENE_OIN_PATH) as f:
            oin_a = f.read().strip()

        self.assertIn("{0>}", oin_a, "Fixture golden must contain a winding marker to flip")

        # Construct a second OIN string with ring 0's winding flipped.
        oin_b = oin_a.replace("{0>}", "{0<}", 1)
        self.assertNotEqual(oin_a, oin_b, "Flip must actually change the input string")

        try:
            reencoded_a = _generate_and_reencode(oin_a)
            reencoded_b = _generate_and_reencode(oin_b)
        except Exception as exc:  # noqa: BLE001 - diagnostic capture
            self.fail(f"generation crashed: {exc!r}")
            return

        winding_re = re.compile(r"\{0[><]\}")
        winding_a = winding_re.search(reencoded_a)
        winding_b = winding_re.search(reencoded_b)

        # Either acceptable diagnostic outcome:
        #  (a) the two flipped inputs re-encode to DIFFERENT winding markers, or
        #  (b) winding markers survive the round trip at all (today: neither
        #      re-encoded string even contains a winding marker on ring 0).
        markers_present = winding_a is not None and winding_b is not None
        markers_differ = markers_present and winding_a.group(0) != winding_b.group(0)

        self.assertTrue(
            markers_present and markers_differ,
            "Winding markers did not survive/distinguish the round trip:\n"
            f"  OIN(A) in:  {oin_a}\n"
            f"  OIN(B) in:  {oin_b}\n"
            f"  OIN(A) out: {reencoded_a}\n"
            f"  OIN(B) out: {reencoded_b}",
        )


class TestHapticFaceCorrectionPhase3(unittest.TestCase):
    """Stereo Phase 3 (Haptic Face Correction) acceptance tests.

    Blessed against the exact rdkit/numpy/scipy versions recorded in the
    module-level _BLESSED_* constants (the entire diagnostic rests on ETKDG
    seed-42 embedding determinism).
    """

    def test_version_pins_match_blessed_golden_environment(self):
        """NFR (SuperPRD Stereo Phase 3, §5.4): record the exact rdkit/numpy/
        scipy versions the seed-42 golden was blessed against, since
        pyproject.toml only pins lower bounds (repo convention). A future
        embedding-handedness change in any of these libraries should surface
        here as a loud, attributable failure rather than a silently wrong
        golden.
        """
        self.assertEqual(rdkit.__version__, _BLESSED_RDKIT_VERSION)
        self.assertEqual(np.__version__, _BLESSED_NUMPY_VERSION)
        self.assertEqual(scipy.__version__, _BLESSED_SCIPY_VERSION)

    @unittest.expectedFailure
    def test_haptic_face_golden_match(self):
        """Test 1 (Deterministic — golden, US-001.3).

        Ferrocene-halide-face.xyz's OIN(1) -> generate -> re-encode should
        reproduce the pinned §5.1 hand-verified golden string exactly.

        Diagnostic result (2026-07-03, post-implementation): the correction
        mechanism itself is verified correct (per-ring winding sense survives
        -- see test_haptic_face_per_ring_flip_inverts_only_that_ring and
        test_haptic_face_two_branch_coverage below, both hard passes), but
        the exact STRING does not match byte-for-byte. Root cause (R2,
        "re-encoder heading-atom instability", anticipated by the SuperPRD):
        the XYZ->OIN re-encoder (i) is free to choose a different specific
        ring atom as the visible heading/marker atom each round trip (here it
        moves between rounds -- confirmed by direct inspection, not a
        one-off), and (ii) is free to list the two (chemically distinct)
        rings in either fragment order. Neither changes the physical winding
        sense or the ring's substituent identity (verified independently by
        the per-ring, content-anchored tests below) -- it is a pre-existing
        re-encoder canonicalization instability, not a generation-side
        correction bug. Kept as a hard `expectedFailure` (never silently
        downgraded) per the MiniPRD's R2 escape hatch, so a future re-encoder
        stability fix flips this to a real pass automatically.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            golden = f.read().strip()

        oin1 = XYZToSMILES().convert(_HALIDE_FACE_XYZ)
        self.assertEqual(oin1, golden, "XYZ->OIN encoding of the fixture must match the golden")

        reencoded = _generate_and_reencode(oin1)
        self.assertEqual(
            golden,
            reencoded,
            f"Generated+re-encoded OIN does not match the pinned golden:\n"
            f"  golden:    {golden}\n"
            f"  reencoded: {reencoded}",
        )

    @unittest.skip(
        "R2 documented fallback (SuperPRD Stereo Phase 3, Risk R2) -- never a "
        "silent auto-downgrade of test_haptic_face_golden_match. This is a "
        "WEAKER, geometric 'the strings differ, but only by a heading-atom / "
        "fragment-order relabeling, never by an actual reversed halide "
        "sequence' check. It is intentionally NOT wired in to replace the "
        "primary exact-match assertion, because a genuine reflection bug "
        "would ALSO present as 'a reversed halide sequence' and must never be "
        "masked by this weaker check. Kept skipped (not deleted) so the "
        "escape hatch is visible and its logic is available, but the exact "
        "assertion above remains the load-bearing one."
    )
    def test_haptic_face_r2_geometric_fallback_never_auto_substituted(self):
        """R2 fallback: verify the golden mismatch is a relabeling, not a
        reversed halide sequence (which would indicate a reflection bug)."""
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            golden = f.read().strip()
        oin1 = XYZToSMILES().convert(_HALIDE_FACE_XYZ)
        reencoded = _generate_and_reencode(oin1)

        # Ring identity via content signature (marker-stripped), NOT position:
        # if the *set* of ring signatures is the same and each ring's own
        # winding character is unchanged, the mismatch above is purely a
        # heading-atom/fragment-order relabeling.
        golden_by_sig = _ring_winding_by_signature(golden)
        reencoded_by_sig = _ring_winding_by_signature(reencoded)
        self.assertEqual(set(golden_by_sig), set(reencoded_by_sig))
        for sig, char in golden_by_sig.items():
            self.assertEqual(
                char,
                reencoded_by_sig[sig],
                f"Ring {sig!r} winding sense changed -- this would indicate a "
                "real reflection/correction bug, not mere relabeling.",
            )

    def test_haptic_face_per_ring_flip_inverts_only_that_ring(self):
        """Test 2 (Deterministic — flip, US-001.2).

        Acceptance test on Ferrocene-halide-face.xyz: per-ring regex
        round-trip anchored to ring CONTENT (substituent signature), not a
        specific heading atom or fragment position (robust to R2). Flipping
        one ring's input winding must invert only that ring's output
        character; the untouched ring's character must be unchanged.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin_a = f.read().strip()

        self.assertIn("{0<}", oin_a, "Fixture golden must contain ring 0's winding marker to flip")
        oin_b = oin_a.replace("{0<}", "{0>}", 1)
        self.assertNotEqual(oin_a, oin_b, "Flip must actually change the input string")

        reencoded_a = _generate_and_reencode(oin_a)
        reencoded_b = _generate_and_reencode(oin_b)

        windings_a = _ring_winding_by_signature(reencoded_a)
        windings_b = _ring_winding_by_signature(reencoded_b)

        self.assertEqual(
            set(windings_a),
            set(windings_b),
            "Same two rings (by content signature) must appear in both outputs:\n"
            f"  A: {reencoded_a}\n  B: {reencoded_b}",
        )
        self.assertTrue(windings_a, "Expected at least one winding marker to survive generation")

        changed = [sig for sig in windings_a if windings_a[sig] != windings_b[sig]]
        unchanged = [sig for sig in windings_a if windings_a[sig] == windings_b[sig]]

        self.assertEqual(
            len(changed),
            1,
            f"Expected exactly one ring's winding to flip:\n  A: {windings_a}\n  B: {windings_b}",
        )
        self.assertEqual(
            len(unchanged),
            1,
            f"Expected exactly one ring's winding to stay put:\n  A: {windings_a}\n  B: {windings_b}",
        )

    def test_haptic_face_two_branch_coverage(self):
        """Test 3 (Deterministic — two branch, US-002.3 / R3).

        The instrumented decision signal (`GeneratedStructure.haptic_face_decisions`)
        must show both branches exercised on the golden fixture: one ring's
        correction already matched (`skipped`) and the other's did not
        (`fired`) -- verified via the structured signal, not by flipping the
        input marker.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin1 = f.read().strip()

        structure = OIN3DGenerator().generate(oin1)
        statuses = sorted(d["status"] for d in structure.haptic_face_decisions)

        self.assertEqual(
            statuses,
            ["fired", "skipped"],
            f"Expected one fired + one skipped ring, got: {structure.haptic_face_decisions}",
        )

    def test_haptic_face_correction_is_proper_rotation(self):
        """Test 4 (Deterministic — det +1, US-006.1).

        The correction rotation built by `_in_plane_correction_axis` +
        `_proper_180_rotation` must be a proper (det +1) rotation. Exercised
        both via the internal helper directly and via the fired ring in the
        golden fixture's generation (no AssertionError raised).
        """
        from oinsmiles.generation.molassembler_adapter import (
            _in_plane_correction_axis,
            _proper_180_rotation,
        )

        binding_pos = np.array(
            [
                [1.0, 0.0, 0.3],
                [0.3, 0.95, 0.1],
                [-0.8, 0.6, -0.1],
                [-0.8, -0.6, 0.05],
                [0.3, -0.95, -0.2],
            ]
        )
        centroid = binding_pos.mean(axis=0)
        axis_unit = np.array([0.0, 0.0, 1.0])
        rot_axis = _in_plane_correction_axis(binding_pos, centroid, axis_unit)
        self.assertIsNotNone(rot_axis)

        rot = _proper_180_rotation(rot_axis)  # raises AssertionError if det != +1
        det = np.linalg.det(rot.as_matrix())
        self.assertAlmostEqual(det, 1.0, places=6)

        # End-to-end: generating the golden fixture (which fires a real
        # correction on ring 0) must not raise the internal det assertion.
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin1 = f.read().strip()
        structure = OIN3DGenerator().generate(oin1)
        fired = [d for d in structure.haptic_face_decisions if d["status"] == "fired"]
        self.assertTrue(fired, "Expected the golden fixture to exercise the fired branch")

    def test_haptic_face_chirality_witness_cip_invariant(self):
        """Test 5 (Deterministic — CIP invariance, US-006.2).

        Chirality-witness fixture: one eta ring (asymmetric, substituted)
        bearing a pendant stereocenter (-CHFCl), plus a second plain/symmetric
        ferrocene-like ring. The pendant stereocenter's CIP code must be
        invariant whether the ring-0 correction fires or is skipped.
        """
        with open(_CHIRALITY_WITNESS_OIN_PATH) as f:
            oin_a = f.read().strip()
        self.assertIn("{0<}", oin_a)
        oin_b = oin_a.replace("{0<}", "{0>}", 1)
        self.assertNotEqual(oin_a, oin_b)

        cips = {}
        for label, oin in (("A", oin_a), ("B", oin_b)):
            structure = OIN3DGenerator().generate(oin)
            decisions = {d["fragment_idx"]: d["status"] for d in structure.haptic_face_decisions}
            mol = structure.mol
            self.assertIsNotNone(mol, f"variant {label}: expected a bonded mol")
            mol.UpdatePropertyCache(strict=False)
            Chem.AssignStereochemistryFrom3D(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

            witness = None
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue
                heavy_nbr_syms = sorted(
                    n.GetSymbol() for n in atom.GetNeighbors() if n.GetAtomicNum() != 1
                )
                n_h = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
                n_h += atom.GetTotalNumHs(includeNeighbors=False)
                if heavy_nbr_syms == ["C", "Cl", "F"] and n_h == 1:
                    witness = atom
                    break
            self.assertIsNotNone(witness, f"variant {label}: could not locate the witness atom")
            cips[label] = witness.GetPropsAsDict().get("_CIPCode")
            self.assertIsNotNone(cips[label], f"variant {label}: witness atom has no CIP code")
            # Sanity: this fixture is designed so ring 0's correction status
            # differs between A (mismatched -> fired) and B (already matches
            # -> skipped); both branches are covered by generate_and_reencode
            # elsewhere, so we only assert the decisions dict is non-empty here.
            self.assertIn(1, decisions)

        self.assertEqual(
            cips["A"],
            cips["B"],
            f"Witness CIP code changed across the correction: A={cips['A']} B={cips['B']}",
        )

    def test_haptic_face_idempotency(self):
        """Test 6 (Deterministic — idempotency).

        generate -> encode -> generate -> encode must converge: each ring's
        winding sense (by content signature, robust to R2 heading-atom
        drift) must be identical after the second round trip.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin1 = f.read().strip()

        reencoded_1 = _generate_and_reencode(oin1)
        reencoded_2 = _generate_and_reencode(reencoded_1)

        windings_1 = _ring_winding_by_signature(reencoded_1)
        windings_2 = _ring_winding_by_signature(reencoded_2)

        self.assertEqual(
            windings_1,
            windings_2,
            "Per-ring winding sense did not converge after a second "
            f"round trip:\n  round 1: {reencoded_1}\n  round 2: {reencoded_2}",
        )

    def test_haptic_face_symmetric_ring_is_inert_noop(self):
        """US-003.3: a detected-symmetric ring is an identity no-op at
        runtime (never a wasted rotation), independent of the marker present.
        """
        with open(_FERROCENE_OIN_PATH) as f:
            oin1 = f.read().strip()
        structure = OIN3DGenerator().generate(oin1)
        self.assertTrue(structure.haptic_face_decisions)
        for d in structure.haptic_face_decisions:
            self.assertTrue(d["symmetric"], d)
            self.assertEqual(d["status"], "no-op", d)

    def test_haptic_face_zero_marker_ring_is_noop(self):
        """Zero-marker eta ring (legacy/hand-authored OIN, no winding
        marker at all): correction is skipped and recorded as `no-op`.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin1 = f.read().strip()
        # Strip both rings' winding markers entirely (zero-marker case).
        oin_no_marker = oin1.replace("{0<}", "{0}").replace("{1<}", "{1}")
        self.assertNotIn("<", oin_no_marker)
        self.assertNotIn(">", oin_no_marker)

        structure = OIN3DGenerator().generate(oin_no_marker)
        self.assertTrue(structure.haptic_face_decisions)
        for d in structure.haptic_face_decisions:
            self.assertIsNone(d["target"])
            self.assertEqual(d["status"], "no-op", d)

    def test_haptic_face_multi_marker_same_slot_raises(self):
        """Multi-marker same-slot (`c{0>}…c{0<}`) is a canonical-form
        violation: raise ValueError, never silently pick a winner.
        """
        with open(_HALIDE_FACE_GOLDEN_PATH) as f:
            oin1 = f.read().strip()
        # Ring 0 currently carries its single marker on the heading O atom
        # ("Oc{0<}1..."); add a second, conflicting marker on another atom
        # in the same ring (slot 0).
        oin_multi_marker = oin1.replace("c{0}(Cl)", "c{0>}(Cl)", 1)
        self.assertNotEqual(oin1, oin_multi_marker)

        with self.assertRaises(ValueError):
            OIN3DGenerator().generate(oin_multi_marker)

    def test_haptic_face_bridged_ansa_conflict_no_regression(self):
        """US-004: bridged ansa-metallocene, conflict path.

        `_ANSA_OIN_CONFLICT` forces ring 0 to want a flip while ring 1 (no
        marker) does not -- a genuine US-004 conflict. Per the negative
        constraint (never an independent per-ring correction inside a
        bridged fragment), placement must be left EXACTLY unchanged: byte
        identical to the committed "before" baseline (captured from the
        current, unmodified `main` prior to this MiniPRD's changes).
        """
        with open(_ANSA_BASELINE_XYZ_PATH) as f:
            baseline_xyz = f.read()

        structure = OIN3DGenerator().generate(_ANSA_OIN_CONFLICT)
        self.assertEqual(structure.xyz, baseline_xyz)

        statuses = {d["status"] for d in structure.haptic_face_decisions}
        self.assertEqual(statuses, {"conflict"}, structure.haptic_face_decisions)

    def test_haptic_face_bridged_ansa_natural_no_regression(self):
        """US-004 companion: the *natural* (unflipped) ansa-metallocene OIN,
        where ring 0 already matches its target (`skipped`) and ring 1 has no
        marker at all. Placement must still equal the committed baseline.
        """
        with open(_ANSA_BASELINE_XYZ_PATH) as f:
            baseline_xyz = f.read()

        structure = OIN3DGenerator().generate(_ANSA_OIN_NATURAL)
        self.assertEqual(structure.xyz, baseline_xyz)
        for d in structure.haptic_face_decisions:
            self.assertEqual(d["status"], "skipped", d)

    def test_haptic_face_bridged_ansa_coherent_fire(self):
        """US-004.1: bridged ansa-metallocene, coherent case.

        `_ANSA_OIN_COHERENT_FIRE` forces BOTH rings to mismatch their target
        winding (both want a flip) -- the coherent (agreeing) branch, as
        opposed to the conflict branch covered above. A single whole-fragment
        correction must fire for both rings (never independently), and must
        not crash or clash (placement still passes the inter-fragment
        collision check inside `_template_generate`).
        """
        structure = OIN3DGenerator().generate(_ANSA_OIN_COHERENT_FIRE)
        statuses = {d["status"] for d in structure.haptic_face_decisions}
        self.assertEqual(statuses, {"fired"}, structure.haptic_face_decisions)
        self.assertEqual(len(structure.haptic_face_decisions), 2)
        # A coherent correction must differ from the (unflipped) "before"
        # baseline -- it actually moved the ligand fragment.
        with open(_ANSA_BASELINE_XYZ_PATH) as f:
            baseline_xyz = f.read()
        self.assertNotEqual(structure.xyz, baseline_xyz)

    def test_haptic_face_emission_assertion_guards_template_gating_hole(self):
        """US-008 guard: an in-scope eta geometry whose `winding_by_slot`
        records a marker, but whose surviving `OINVector`s do not carry it
        (the `oin_parser.py:495` template-gating hole), must fail loudly
        rather than silently behave as a legitimate zero-marker ring.
        """
        from oinsmiles.generation.molassembler_adapter import _template_generate
        from oinsmiles.generation.oin_parser import TEMPLATES, OINVector, ParsedOIN

        geo_code = "LIN"
        template = TEMPLATES[geo_code]
        slot0_vec = tuple(float(x) for x in template[0])

        vectors = [
            OINVector(
                atom_idx=-1, vector=slot0_vec, fragment_idx=1, atom_in_fragment_idx=i, winding=None
            )
            for i in range(5)
        ]
        parsed = ParsedOIN(
            smiles="[Fe].[cH]1[cH][cH][cH][cH]1",
            fragments=["[Fe]", "[cH]1[cH][cH][cH][cH]1"],
            metal_fragment_idx=0,
            vectors=vectors,
            original_oin="[Fe_LIN].[cH]1[cH][cH][cH][cH]1",
            geo_code=geo_code,
            # Claims a marker existed for slot 0, but none of `vectors` above
            # carries it -- simulates the oin_parser.py hole.
            winding_by_slot={0: "<"},
        )

        with self.assertRaises(AssertionError):
            _template_generate(parsed)


if __name__ == "__main__":
    unittest.main()
