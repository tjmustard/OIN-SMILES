"""The guard that stops OIN_STABLE_STEREO from "fixing" the bug by deleting the answer.

A descriptor that is stable under atom renumbering *because it is constant* is worse than
the instability it replaces: the string stops distinguishing enantiomers and nothing reports
it. That is not hypothetical -- the Y2 axial wave shipped exactly that (a token sorted by
sign, which made it reflection-invariant), and every guard written against the easy
single-axis fixture passed. Only a mirror check caught it.

So stability is only half the requirement. These tests assert BOTH halves:

    stable   the same structure, atoms renumbered, encodes byte-identically
    faithful its z-mirror encodes DIFFERENTLY, and specifically to the same string with
             every chiral tag inverted -- a textbook enantiomer pair

The second is the one that has teeth. Counting tags is NOT sufficient to check it: a
molecule with three ``@@`` and three ``@`` mirrors to three ``@`` and three ``@@``, which has
identical counts. These tests compare the whole string under an ``@``<->``@@`` swap instead.
"""

import os
import random
import re
import tempfile
import unittest

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "fixtures")

# Both carry stereocentres the encoder captures, and both were in the measured set of 29
# molecules whose absolute stereochemistry flipped under pure atom renumbering.
STEREO_FIXTURES = ["EJUJUP_comp_0.xyz", "OCUGIC_comp_0.xyz"]


def _swap_tags(s):
    """Exchange ``@`` and ``@@`` without double-substituting."""
    return re.sub(r"@@|@(?!@)", lambda m: "@" if m.group(0) == "@@" else "@@", s)


def _read(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].split()[0])
    rows = [ln.split() for ln in lines[2 : 2 + n]]
    return [(r[0], float(r[1]), float(r[2]), float(r[3])) for r in rows], lines[1]


def _write(rows, comment):
    path = os.path.join(tempfile.mkdtemp(), "s.xyz")
    with open(path, "w") as f:
        f.write(f"{len(rows)}\n{comment if comment.endswith(chr(10)) else comment + chr(10)}")
        for sym, x, y, z in rows:
            f.write(f"{sym:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")
    return path


def _available():
    return [f for f in STEREO_FIXTURES if os.path.exists(os.path.join(FIXTURES, f))]


@unittest.skipUnless(_available(), "stereo-flip regression fixtures not present")
class TestStableStereoIsBothStableAndFaithful(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from oinsmiles import XYZToSMILES

        cls.conv = XYZToSMILES()
        cls._prior = os.environ.get("OIN_STABLE_STEREO")
        os.environ["OIN_STABLE_STEREO"] = "1"

    @classmethod
    def tearDownClass(cls):
        if cls._prior is None:
            os.environ.pop("OIN_STABLE_STEREO", None)
        else:
            os.environ["OIN_STABLE_STEREO"] = cls._prior

    def test_renumbering_is_byte_stable(self):
        """Half one: presentation must not change the string."""
        for fx in _available():
            path = os.path.join(FIXTURES, fx)
            rows, comment = _read(path)
            base = self.conv.convert(path)
            self.assertIsNotNone(base, fx)
            for seed in (1, 7, 13):
                shuffled = rows[:]
                random.Random(seed).shuffle(shuffled)
                got = self.conv.convert(_write(shuffled, comment))
                self.assertEqual(base, got, f"{fx}: renumbering (seed {seed}) changed the string")

    def test_mirror_does_not_collapse(self):
        """Half two, the one with teeth: the descriptor must not be constant.

        The assertion is deliberately about WHERE the difference lives, not about a specific
        tag flip, because this notation carries stereochemistry in two places: sp3 chiral
        tags (``@``/``@@``) and the eta winding character (``>``/``<``, the coordinated ring
        face). Mirroring an eta complex can flip the winding rather than a chiral tag --
        measured on ``EJUJUP_comp_0``, a Cr arene, whose mirror differs from the original
        only in ``{0>}`` -> ``{0<}``.

        So: the mirror must differ (no collapse), and once every stereo-bearing token is
        neutralized the two must become identical -- i.e. the difference is confined to
        stereochemistry and did not smuggle in a constitution change.
        """
        neutralize = lambda s: re.sub(r"@@|@(?!@)", "", s).replace(">", "").replace("<", "")  # noqa: E731
        for fx in _available():
            path = os.path.join(FIXTURES, fx)
            rows, comment = _read(path)
            base = self.conv.convert(path)
            mirrored = [(s, x, y, -z) for (s, x, y, z) in rows]
            mirror = self.conv.convert(_write(mirrored, comment))

            self.assertIsNotNone(mirror, fx)
            self.assertTrue(
                ("@" in base) or (">" in base) or ("<" in base),
                f"{fx}: fixture carries no stereo token at all, so it cannot guard anything",
            )
            self.assertNotEqual(
                base, mirror, f"{fx}: mirror COLLAPSED onto the original -- stereo was destroyed"
            )
            self.assertEqual(
                neutralize(base),
                neutralize(mirror),
                f"{fx}: mirror differs OUTSIDE the stereo tokens, so mirroring changed "
                "constitution rather than only configuration",
            )

    def test_pure_sp3_case_inverts_every_chiral_tag(self):
        """Stricter form, for a fixture whose stereo is sp3 tags only.

        Here the mirror must be exactly the base with ``@``<->``@@`` exchanged. Tag COUNTS
        are not sufficient to check this -- three ``@@`` plus three ``@`` mirrors to three
        ``@`` plus three ``@@``, identical counts -- so the whole string is compared under
        the swap.
        """
        fx = "OCUGIC_comp_0.xyz"
        if not os.path.exists(os.path.join(FIXTURES, fx)):
            self.skipTest(f"{fx} not present")
        path = os.path.join(FIXTURES, fx)
        rows, comment = _read(path)
        base = self.conv.convert(path)
        if ">" in base or "<" in base:
            self.skipTest(f"{fx} carries eta winding; covered by test_mirror_does_not_collapse")
        mirror = self.conv.convert(_write([(s, x, y, -z) for (s, x, y, z) in rows], comment))
        self.assertEqual(_swap_tags(base), mirror, f"{fx}: not encoded as an enantiomer pair")


if __name__ == "__main__":
    unittest.main()
