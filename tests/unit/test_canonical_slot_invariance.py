"""End-to-end canonicality guards for ``OIN_CANONICAL_SLOTS`` (v0.4.5 Lane 2).

``test_canonical_slots.py`` tests the machinery on hand-written strings. This file tests the
**encoder**: it re-presents a real fixture -- permuted atom lines, a random proper rotation,
or both -- and asserts one byte-identical OIN string, because the molecular graph never
changed so the correct answer is byte-identical rather than merely similar.

Two things are deliberately asserted together:

* with the lever ON, drift goes away;
* with the lever OFF, output is byte-identical to the pinned goldens
  (``test_regression_stability.py`` owns the full golden set; this file re-checks the two
  that exercise the seam most directly).

The rotations are forced **proper** (det = +1). An improper operation mirrors the structure,
which legitimately changes a chiral molecule's encoding, so folding it in would be testing
the wrong invariance -- the same reason ``derive_rotation_group`` filters on ``det > 0``.
"""

import os
import random
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES  # noqa: E402

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))

#: Small, fast fixtures that between them cover the seams: a homoleptic square plane
#: (slot relabeling only), a chelate (donor pair inside one fragment), a trigonal
#: bipyramid (a geometry whose group is not the octahedron's), and a mixed-donor square
#: plane whose four ligands are all distinct (no automorphism to hide behind).
_CASES = ("CisPlatin.xyz", "Cis-PtCl2(en).xyz", "FeCO5.xyz", "PtMeNH3ClBr-Cis.xyz")


def _read_xyz(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.asarray(coords, dtype=float)


def _write_xyz(path, syms, coords):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n\n")
        for sym, xyz in zip(syms, coords):
            fh.write(f"{sym:<3} {xyz[0]:>16.10f} {xyz[1]:>16.10f} {xyz[2]:>16.10f}\n")


def _proper_rotation(rng):
    """A random rotation with det = +1 -- never a reflection of the result."""
    a = np.asarray([[rng.gauss(0, 1) for _ in range(3)] for _ in range(3)], dtype=float)
    q, r = np.linalg.qr(a)
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


class _PresentationMixin:
    @classmethod
    def setUpClass(cls):
        cls.conv = XYZToSMILES()
        cls.tmpdir = tempfile.mkdtemp(prefix="oin-lane2-")

    def _presentations(self, name, trials=3, seed=11):
        """Yield ``(label, path)`` for renumbered / rotated / both re-presentations."""
        syms, coords = _read_xyz(os.path.join(_FIXTURES, name))
        rng = random.Random(seed)
        for trial in range(trials):
            for mode in ("renumber", "rotate", "both"):
                order = list(range(len(syms)))
                out = coords
                if mode in ("renumber", "both"):
                    rng.shuffle(order)
                    out = coords[order]
                if mode in ("rotate", "both"):
                    out = out @ _proper_rotation(rng).T
                path = os.path.join(self.tmpdir, f"{name}.{mode}.{trial}.xyz")
                _write_xyz(path, [syms[i] for i in order], out)
                yield f"{mode}/{trial}", path


class TestCanonicalSlotsAreInvariant(_PresentationMixin, unittest.TestCase):
    """With the lever on, presentation must not reach the string."""

    def test_byte_identical_under_renumbering_and_rotation(self):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_SLOTS": "1"}):
            for name in _CASES:
                base = self.conv.convert(os.path.join(_FIXTURES, name))
                for label, path in self._presentations(name):
                    self.assertEqual(
                        self.conv.convert(path),
                        base,
                        f"{name} drifted under {label}",
                    )

    def test_renumbering_alone_is_the_hard_case(self):
        """Documents why the fix is needed: rotation was already fine, numbering was not.

        Asserted as a *property of the lever*, not of the baseline -- the baseline's drift
        is measured in ``docs/RENUMBERING_INSTABILITY_v0.4.5.md`` and is molecule-dependent,
        so pinning it here would be pinning a number that is allowed to improve.
        """
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_SLOTS": "1"}):
            name = "PtMeNH3ClBr-Cis.xyz"
            base = self.conv.convert(os.path.join(_FIXTURES, name))
            seen = {
                self.conv.convert(path)
                for label, path in self._presentations(name, trials=4)
                if label.startswith("renumber")
            }
            self.assertEqual(seen, {base})


class TestLeverOffIsByteIdentical(unittest.TestCase):
    """Unset ``OIN_CANONICAL_SLOTS`` must reproduce the pre-Lane-2 bytes exactly.

    The post-pass, the ``_align_to_pai`` neutralization and the fragment tie-break are all
    behind the one lever; the rotation-group unification is not, because it is provably a
    no-op on every geometry except PBP (see ``test_canonical_slots.py``).
    """

    def test_goldens(self):
        # "0", not absent: the lever is default-ON since v0.4.5, so deleting the key would
        # exercise the ON path and then assert it produced the OFF bytes.
        env = dict(os.environ, OIN_CANONICAL_SLOTS="0")
        with mock.patch.dict(os.environ, env, clear=True):
            conv = XYZToSMILES()
            self.assertEqual(
                conv.convert(os.path.join(_FIXTURES, "CisPlatin.xyz")),
                "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
            )
            self.assertEqual(
                conv.convert(os.path.join(_FIXTURES, "fac-Ir(ppy)3.xyz")),
                "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1",
            )


class TestLeverOnGoldens(unittest.TestCase):
    """Pins WHAT the lever changes, and proves each change is a relabeling, not a re-isomer.

    Measured with ``OIN_CANONICAL_SLOTS`` alone over the six pinned fixtures: three move
    (CisPlatin, TransPlatin, fac-Ir(ppy)3) and three are already canonical (mer-Ir(ppy)3,
    Ferrocene, Cis-PtCl2(en)). Pinning the moved strings turns the promotion A/B into a diff
    against a committed expectation rather than a judgement call; pinning the *unmoved* ones
    catches gratuitous churn.

    Cisplatin is the clearest case. The old labeling put the chlorides on slots 0,1 only
    because the fragment sort ranks by descending mass. The lex-min runs over the vertex
    COLOURS and ``"N" < "[Cl]"`` bytewise, so the canonical labeling puts the amines on the
    lowest vertices. Cis is still cis: the chlorides land on 2,3, which are adjacent.

    The second assertion in each case is the one that matters -- the canonical round-trip key
    is UNCHANGED, i.e. the lever relabeled the isomer without turning it into a different
    one. A canonicalization that merged two isomers would fail here even with the string
    assertion passing.
    """

    _EXPECT = {
        "CisPlatin.xyz": "[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}",
        "TransPlatin.xyz": "[Pt_SPL].N{0}.[Cl]{1}.N{2}.[Cl]{3}",
        "fac-Ir(ppy)3.xyz": (
            "[Ir_OCT].c{0}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{1}1.c{4}1ccccc1-c1ccccn{3}1"
        ),
    }
    _UNMOVED = ("mer-Ir(ppy)3.xyz", "Ferrocene.xyz", "Cis-PtCl2(en).xyz")

    @staticmethod
    def _encode(name, lever_on):
        # Both directions explicit -- see TestLeverOffIsByteIdentical.test_goldens.
        env = dict(os.environ, OIN_CANONICAL_SLOTS="1" if lever_on else "0")
        with mock.patch.dict(os.environ, env, clear=True):
            return XYZToSMILES().convert(os.path.join(_FIXTURES, name))

    def test_moved_goldens_are_relabelings_not_new_isomers(self):
        from oinsmiles.oin.compare import canonical_roundtrip_key

        for name, expected in self._EXPECT.items():
            with self.subTest(name):
                off, on = self._encode(name, False), self._encode(name, True)
                self.assertNotEqual(off, on, "expected this fixture's labeling to move")
                self.assertEqual(on, expected, "the canonical labeling itself changed")
                self.assertEqual(
                    canonical_roundtrip_key(off),
                    canonical_roundtrip_key(on),
                    "the lever changed the ISOMER, not just the labels -- over-folding",
                )

    def test_already_canonical_fixtures_do_not_churn(self):
        for name in self._UNMOVED:
            with self.subTest(name):
                self.assertEqual(self._encode(name, False), self._encode(name, True))


class TestNoOverFolding(unittest.TestCase):
    """The whole risk of this lane: a canonicalization that merges two real isomers.

    Drift is visible; a silent merge is not. ``test_facmer_key.py`` and
    ``tests/integration/test_isomer_divergence.py`` are the primary guards; this adds the
    end-to-end check that the *emitted string* still separates them with the lever ON.
    """

    def _encode(self, name):
        with mock.patch.dict(os.environ, {"OIN_CANONICAL_SLOTS": "1"}):
            return XYZToSMILES().convert(os.path.join(_FIXTURES, name))

    def test_fac_and_mer_stay_distinct(self):
        self.assertNotEqual(
            self._encode("fac-Ir(ppy)3.xyz"),
            self._encode("mer-Ir(ppy)3.xyz"),
        )

    def test_cis_and_trans_stay_distinct(self):
        self.assertNotEqual(
            self._encode("PtMeNH3ClBr-Cis.xyz"),
            self._encode("PtMeNH3ClBr-Trans.xyz"),
        )

    def test_cisplatin_and_transplatin_stay_distinct(self):
        self.assertNotEqual(self._encode("CisPlatin.xyz"), self._encode("TransPlatin.xyz"))


if __name__ == "__main__":
    unittest.main()
