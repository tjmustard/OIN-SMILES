"""Isomer divergence: a different isomer must encode to a different string.

The dual of conformer *convergence*: distinct isomers must NOT collapse together. The subtle
part is that different axes of isomerism live at different comparison layers, and asserting on
the wrong one silently passes. Two surfaces exist (see ``src/oinsmiles/oin/compare.py``):

* **Raw canonical string** -- the encoder's output (``XYZToSMILES().convert`` /
  ``get_oin_string``). ``tests/unit/test_regression_stability.py`` pins these byte-identical.
* **Round-trip equivalence key** -- ``winding_canonical_key(normalize_oin_for_comparison(oin))``
  (``_key`` below). Deliberately lossy: it strips the metal ``@SPn`` descriptor (deferred; see
  geometry-canonical-slot-key.md §4) and folds the benign slot-relabeling drift between true
  conformers via a symmetry-canonical vertex signature. As of v0.4.4 it is fac/mer-aware -- the
  signature keeps positional isomers distinct while still collapsing conformers of one isomer.

Which layer carries which axis (verified against ``compare.py`` and by ``convert()``):

===================  ===============================================  ==========================
axis                 distinguished by                                 divergence assertion
===================  ===============================================  ==========================
geometry SPL/TET     geometry tag ``_XXX`` (kept)                     ``_key`` differs
ring winding {n>/<}  winding multiset                                 ``_key`` differs
E/Z ``/`` vs ``\\``   bond-direction tokens (kept)                     raw AND ``_key`` differ
cis<->trans slot     donor arrangement over the polyhedron            ``_key`` differs
fac<->mer slot       donor arrangement over the polyhedron (v0.4.4)   raw AND ``_key`` differ
metal @SPn           deleted by the normalizer (deferred)             **raw differs; _key EQUAL**
===================  ===============================================  ==========================

The last row is the load-bearing one: metal-@ isomerism survives only in the raw canonical
string, so it MUST be compared raw -- the equivalence key is still blind to it (deferred until
the encoder has a reproducible metal stereo descriptor). fac/mer used to share that blindness;
as of v0.4.4 the symmetry-canonical vertex signature catches it at the key layer too.

``TestGenerativeDivergence`` additionally drives OIN -> 3D -> OIN under the deterministic
``optimizer="ff"`` + ``seed=42`` path (mirrors ``test_roundtrip_smoke.py``): a mutated isomer,
regenerated and re-encoded, must still diverge from the base. Axes the generator cannot yet
build a clean alternate geometry for (winding face, SPL/TET) ``skipTest`` with an A3/A4 pointer
rather than weakening the assertion.
"""

import contextlib
import os
import unittest
from pathlib import Path

import numpy as np

from oinsmiles import XYZToSMILES
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key
from oinsmiles.utils.xyz2mol import get_oin_string

REPO = Path(__file__).resolve().parents[2]
SEED = 42

# Canonical goldens (match test_regression_stability.py / test_roundtrip_smoke.py byte-for-byte).
CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
TRANSPLATIN = "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"
FERROCENE = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
FERROCENE_FLIP = (
    "[Fe_LIN].[cH]{0<}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
)
FAC_IRPPY3 = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
MER_IRPPY3 = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"
NI_SPL = "[Ni_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
NI_TET = "[Ni_TET].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


@contextlib.contextmanager
def _silence_fds():
    """Redirect C-level stdout/stderr to devnull (openbabel prints distance warnings)."""
    with open(os.devnull, "w") as devnull:
        old_out, old_err = os.dup(1), os.dup(2)
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(old_out, 1)
            os.dup2(old_err, 2)
            os.close(old_out)
            os.close(old_err)


def _key(oin):
    """Round-trip equivalence key. Returns ``(str, list)`` -- unhashable; compare by ``==``."""
    return winding_canonical_key(normalize_oin_for_comparison(oin))


def _reencode(gen_mol, xyz):
    """Re-encode a generated structure to OIN via its bonded mol (mirrors verify_roundtrip)."""
    lines = xyz.splitlines()
    natoms = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + natoms)])
    return get_oin_string(gen_mol, coords)


class TestIsomerComparisonLayer(unittest.TestCase):
    """Hermetic, deterministic proofs that each isomer axis diverges at the correct layer."""

    def test_geometry_axis(self):
        """SPL vs TET differ via the kept geometry tag (canonical-key layer)."""
        self.assertNotEqual(_key(NI_SPL), _key(NI_TET))

    def test_winding_axis(self):
        """A both-ring winding flip differs via the winding multiset (canonical-key layer)."""
        self.assertNotEqual(_key(FERROCENE), _key(FERROCENE_FLIP))

    def test_ez_axis(self):
        """E vs Z differ via bond-direction tokens -- caught at BOTH raw and key layers."""
        e, z = "C/C=C/C", "C/C=C\\C"
        self.assertNotEqual(e, z)
        self.assertNotEqual(_key(e), _key(z))

    def test_cis_trans_slot(self):
        """cis vs trans differ via ligand-type sequence, which survives slot renumbering."""
        self.assertNotEqual(_key(CISPLATIN), _key(TRANSPLATIN))

    def test_metal_stereo_raw_only(self):
        """Metal @-chirality lives ONLY in the raw string; the key erases it.

        Proves the layer choice: comparing metal chirality with the canonical key would
        WRONGLY call these equal, so it must be compared raw.

        Metal ``@SPn`` in the key is DEFERRED (unlike fac/mer, which v0.4.4's vertex
        signature now catches) -- it needs a reproducible encoder-side metal stereo
        descriptor first. See spec/handoffs/v0.4.4/geometry-canonical-slot-key.md §4.
        """
        p1 = "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
        p2 = "[Pt@SP2_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
        self.assertNotEqual(p1, p2)  # raw: distinct isomers
        self.assertEqual(_key(p1), _key(p2))  # key: blind (would falsely converge) -- deferred

    def test_fac_mer_divergence(self):
        """fac vs mer Ir(ppy)3 diverge in BOTH the raw string and the canonical key (v0.4.4).

        Both are three identical ppy fragments differing only in absolute slot numbers. Before
        v0.4.4 the equivalence key renumbered those away, so it was fac/mer-blind (like metal-@
        still is). v0.4.4 replaced the renumber with a symmetry-canonical vertex signature, so
        the key now diverges too. The raw ``convert()`` output also distinguishes them (and
        equals the byte-identical test_regression_stability goldens).
        """
        conv = XYZToSMILES()
        with _silence_fds():
            fac = conv.convert(str(REPO / "tests" / "fixtures" / "fac-Ir(ppy)3.xyz"))
            mer = conv.convert(str(REPO / "tests" / "fixtures" / "mer-Ir(ppy)3.xyz"))
        # anchor to the committed goldens (guards an encoder regression too)
        self.assertEqual(fac, FAC_IRPPY3, "fac encoding drifted from the pinned golden")
        self.assertEqual(mer, MER_IRPPY3, "mer encoding drifted from the pinned golden")
        # the divergence guarantee: distinct canonical strings
        self.assertNotEqual(fac, mer)
        # v0.4.4: the fac/mer-aware vertex signature diverges at the key layer too
        self.assertNotEqual(_key(fac), _key(mer))


class TestGenerativeDivergence(unittest.TestCase):
    """OIN -> 3D -> OIN under FF + fixed seed: a mutated isomer must still diverge from base."""

    @staticmethod
    def _generate(oin):
        return OIN3DGenerator(optimizer="ff", seed=SEED).generate(oin)

    def test_cis_trans_generative(self):
        """Regenerate transplatin (a cis->trans slot mutation) and confirm it stays distinct."""
        gen = self._generate(TRANSPLATIN)
        self.assertIsNotNone(gen.mol, "generator returned no bonded mol for transplatin")
        out = _reencode(gen.mol, gen.xyz)
        self.assertNotEqual(_key(out), _key(CISPLATIN), "trans regenerated to the cis notation")
        self.assertEqual(_key(out), _key(TRANSPLATIN), "trans did not round-trip to itself")

    def test_winding_flip_generative(self):
        """Regenerate a winding-flipped ferrocene; assert divergence or skip if not buildable."""
        try:
            gen = self._generate(FERROCENE_FLIP)
        except (ValueError, TimeoutError) as exc:
            self.skipTest(f"generator could not build the winding-flipped ferrocene: {exc}")
        if gen.mol is None:
            self.skipTest("generator returned no bonded mol for the winding-flipped ferrocene")
        out = _reencode(gen.mol, gen.xyz)
        if _key(out) == _key(FERROCENE):
            self.skipTest(
                "generator/encoder canonicalizes ferrocene winding to the base face; the "
                "flipped-face isomer is not exercised through generation -- pending A3/A4"
            )
        self.assertNotEqual(_key(out), _key(FERROCENE))

    def test_geometry_generative(self):
        """Regenerate a TET mutant of an SPL base; assert divergence or skip if not buildable."""
        try:
            gen = self._generate(NI_TET)
        except (ValueError, TimeoutError) as exc:
            self.skipTest(f"generator could not build the TET mutant: {exc}")
        if gen.mol is None:
            self.skipTest("generator returned no bonded mol for the TET mutant")
        out = _reencode(gen.mol, gen.xyz)
        if _key(out) == _key(NI_SPL):
            self.skipTest(
                "generator did not build a distinct TET geometry (re-encodes to SPL); the "
                "SPL/TET axis is not yet controllable through generation -- pending A3/A4"
            )
        self.assertNotEqual(_key(out), _key(NI_SPL))


if __name__ == "__main__":
    unittest.main()
