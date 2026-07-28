"""The reflection-parity veto for the donor fold (v0.4.12 Lane 1, ``OIN_FOLD_PARITY_VETO``).

WHY THESE TESTS ARE AT THE COORDINATE LEVEL AND NOT THE STRING LEVEL
====================================================================
``tests/unit/test_canonical_slots.py::TestDonorFoldCollapsesEnantiomers`` pins the defect on
two hand-written OIN strings and **still passes unchanged** after this lever landed. That is
not an oversight, and it is the single most transferable fact in this file:

    Reflection parity is not a property of the emitted string. Two labelings that differ by a
    donor swap are related by a permutation whose PARITY depends on where those donors sit in
    space, and the string does not carry that.

So the fix cannot live where the defect is visible. ``canonicalize_oin_slots`` collapses the
pair exactly as before; the veto sits one level up, in ``get_oin_string``, where the pristine
conformer is still in hand. Every test that could show the fix therefore has to start from
coordinates.

The three fixtures are the molecules v0.4.11 confirmed genuinely chiral with
``tools/injectivity/oracle.py`` **without hitting its automorphism cap**, so the verdict does
not rest on the unreliable ones.
"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from oinsmiles import XYZToSMILES
from oinsmiles.oin import fold_parity
from oinsmiles.oin.levers import held_off, lever_enabled

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "fold_parity"

#: Confirmed chiral, and confirmed collapsed by the bare fold (v0.4.11 Lane 2 §3).
COLLAPSE_FIXTURES = ["BIWDIV_comp_0", "CIHVAT_comp_0", "OJEKET_comp_0"]

FOLD = "OIN_CANONICAL_DONOR_FOLD"
VETO = "OIN_FOLD_PARITY_VETO"


def _read(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms = [ln.split()[0] for ln in lines[2 : 2 + n]]
    xyz = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines[2 : 2 + n]])
    return syms, xyz


def _write(syms, xyz):
    fh = tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False)
    fh.write(f"{len(syms)}\n\n")
    for s, c in zip(syms, xyz):
        fh.write(f"{s:<3} {c[0]:>16.10f} {c[1]:>16.10f} {c[2]:>16.10f}\n")
    fh.close()
    return fh.name


def _encode(path, fold, veto):
    with mock.patch.dict(os.environ, {FOLD: "1" if fold else "0", VETO: "1" if veto else "0"}):
        return XYZToSMILES().convert(str(path))


def _pair(name, fold, veto):
    """``(encoding, encoding-of-its-mirror)`` under one lever configuration."""
    src = FIXTURES / f"{name}.xyz"
    syms, xyz = _read(src)
    mirrored = xyz.copy()
    mirrored[:, 2] *= -1.0
    mpath = _write(syms, mirrored)
    try:
        return _encode(src, fold, veto), _encode(mpath, fold, veto)
    finally:
        os.unlink(mpath)


class TestDefaultOff(unittest.TestCase):
    def test_registered_and_held_off(self):
        self.assertIn(VETO, held_off())

    def test_unset_means_off(self):
        """The trap that cost 23 test failures across two promotions."""
        env = {k: v for k, v in os.environ.items() if k != VETO}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertFalse(lever_enabled(VETO))

    def test_explicit_zero_disables(self):
        with mock.patch.dict(os.environ, {VETO: "0"}):
            self.assertFalse(lever_enabled(VETO))


class TestTheVetoSeparatesConfirmedEnantiomers(unittest.TestCase):
    """The whole point of the lever, on the three oracle-confirmed molecules."""

    def test_shipped_encoder_already_separates_them(self):
        """Control. If this fails the fixture is not a mirror pair and nothing below means anything."""
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                a, b = _pair(name, fold=False, veto=False)
                self.assertNotEqual(a, b, f"{name}: shipped encoder must keep these apart")

    def test_the_bare_fold_collapses_them(self):
        """🔴 The v0.4.11 refutation, reproduced from coordinates rather than from strings."""
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                a, b = _pair(name, fold=True, veto=False)
                self.assertEqual(a, b, f"{name}: the bare fold is expected to collapse this pair")

    def test_the_veto_restores_the_separation(self):
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                a, b = _pair(name, fold=True, veto=True)
                self.assertNotEqual(a, b, f"{name}: the veto must keep these enantiomers apart")

    def test_it_separates_them_BY_VETOING_and_not_by_giving_up(self):
        """🔴 The assertion above passed for the WRONG REASON once, and would again.

        Emitting the rotation-only labeling separates a mirror pair whether the veto reasoned
        about parity or simply had no evidence and declined. The first implementation of this
        module encoded the mirror with the fold INHERITED (on) rather than forced off, which
        made ``s_rot_m`` identical to ``s_fold_m``, disarmed the achiral guard, and tripped the
        self-check on **18 of 18** movers -- and every test in this class still passed.

        So the outcome, not just the string, is the assertion.

        For BIWDIV the fold does not fire on the deposited structure at all -- it fires on its
        MIRROR, moving the mirror's labeling onto the original's. So the assertion is a
        property of the PAIR, not of either member: at least one side must be vetoed on the
        evidence, and no side may reach its answer by declining.
        """
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                outcomes = []
                src = FIXTURES / f"{name}.xyz"
                syms, xyz = _read(src)
                mirrored = xyz.copy()
                mirrored[:, 2] *= -1.0
                mpath = _write(syms, mirrored)
                try:
                    for path in (src, mpath):
                        _encode(path, fold=True, veto=True)
                        outcomes.append(fold_parity.last_outcome())
                finally:
                    os.unlink(mpath)

                self.assertIn(
                    "vetoed_collapse",
                    outcomes,
                    f"{name}: the pair must separate because the veto FIRED, not because the "
                    f"veto gave up. Got {outcomes}.",
                )
                declines = [o for o in outcomes if o and o.startswith("declined")]
                self.assertEqual(
                    declines, [], f"{name}: a decline means there was no instrument, not a verdict"
                )


class TestDefaultPathUntouched(unittest.TestCase):
    def test_veto_alone_changes_nothing(self):
        """With the fold off there is no widened candidate set to police."""
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                src = FIXTURES / f"{name}.xyz"
                self.assertEqual(
                    _encode(src, fold=False, veto=False),
                    _encode(src, fold=False, veto=True),
                    f"{name}: the veto must be inert unless the fold is on",
                )


class TestReconstructionFidelity(unittest.TestCase):
    """``tmc_mol``'s atom order is NOT the coordinate order, and getting that wrong is silent."""

    def test_atom_coord_pairs_honours_origidx(self):
        from oinsmiles.utils.perception_tmc import get_tmc_mol

        src = FIXTURES / "BIWDIV_comp_0.xyz"
        tmc_mol, coords = get_tmc_mol(src, 0, with_stereo=False)
        pairs = fold_parity._atom_coord_pairs(tmc_mol, coords)
        self.assertIsNotNone(pairs)

        # Zipping positionally instead produced [Co_TBP] with invented bond orders in place of
        # [Co_OCT] -- a chemically different molecule that encodes perfectly cleanly. The
        # reconstruction must reproduce the ELEMENT MULTISET of the input file exactly.
        syms, _ = _read(src)
        self.assertEqual(sorted(s for s, _c in pairs), sorted(syms))

    def test_a_mismatched_count_declines_rather_than_guesses(self):
        from oinsmiles.utils.perception_tmc import get_tmc_mol

        tmc_mol, coords = get_tmc_mol(FIXTURES / "BIWDIV_comp_0.xyz", 0, with_stereo=False)
        self.assertIsNone(fold_parity._atom_coord_pairs(tmc_mol, coords[:-1]))


class TestReentrancy(unittest.TestCase):
    """The veto encodes a mirror, and that encode runs through the veto's own call site."""

    def test_guard_is_clear_outside_a_mirror_pass(self):
        self.assertFalse(fold_parity._in_mirror_pass())

    def test_guard_is_set_inside_and_cleared_after(self):
        with fold_parity._mirror_pass():
            self.assertTrue(fold_parity._in_mirror_pass())
        self.assertFalse(fold_parity._in_mirror_pass())

    def test_resolve_is_a_plain_postpass_inside_a_mirror_pass(self):
        """Recursion terminates because the inner encode never re-enters the veto."""
        from oinsmiles.oin.canonical_slots import canonicalize_oin_slots

        oin = "[Co_OCT].O=C(N{0}c1ccccn1)c1cccc(C(=O)N{1}c2ccccn2)n{2}1"
        with mock.patch.dict(os.environ, {FOLD: "1", VETO: "1"}):
            with fold_parity._mirror_pass():
                self.assertEqual(fold_parity.resolve(oin, None, None), canonicalize_oin_slots(oin))

    def test_forced_lever_restores_an_unset_variable_by_deleting_it(self):
        env = {k: v for k, v in os.environ.items() if k != FOLD}
        with mock.patch.dict(os.environ, env, clear=True):
            with fold_parity._forced_lever(FOLD, False):
                self.assertEqual(os.environ[FOLD], "0")
            self.assertNotIn(FOLD, os.environ, "unset and '0' are different states")


class TestPresentationInvariance(unittest.TestCase):
    """The verdict must be a property of the STRUCTURE, not of the input atom numbering."""

    def test_renumbering_the_input_does_not_move_the_answer(self):
        rng = np.random.default_rng(42)
        for name in COLLAPSE_FIXTURES:
            with self.subTest(name):
                src = FIXTURES / f"{name}.xyz"
                syms, xyz = _read(src)
                order = rng.permutation(len(syms))
                shuffled = _write([syms[i] for i in order], xyz[order])
                try:
                    self.assertEqual(
                        _encode(src, fold=True, veto=True),
                        _encode(shuffled, fold=True, veto=True),
                        f"{name}: the veto's verdict must not depend on input numbering",
                    )
                finally:
                    os.unlink(shuffled)


if __name__ == "__main__":
    unittest.main()
