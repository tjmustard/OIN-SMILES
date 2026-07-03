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

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator

_FIXTURES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_CANDIDATES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../candidate_outputs"))

_BDPP_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDPP.xyz")
_BDNN_XYZ = os.path.join(_FIXTURES_DIR, "PdCl2-RR-BDNN.xyz")
_FERROCENE_OIN_PATH = os.path.join(_CANDIDATES_DIR, "ferrocene_oin.txt")


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

    @unittest.expectedFailure
    def test_haptic_face_winding(self):
        """Fixed by Phase 3: eta-ring winding markers ({n>} vs {n<}) must
        survive OIN -> 3D generation -> re-encoding, and flipping one ring's
        winding in the input OIN must produce a distinguishable output.

        Diagnostic result (2026-07-03): winding markers do not survive
        generation at all (SLOT_REGEX only captures the numeric slot rank,
        discarding the direction suffix) — see NOTES.md Log entry.
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


if __name__ == "__main__":
    unittest.main()
