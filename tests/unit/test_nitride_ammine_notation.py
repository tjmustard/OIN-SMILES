"""A terminal (heavy==0) nitride donor must be notationally distinguishable from an
ammine in the OIN string, or the round trip re-protonates it.

Both a bare nitride ``[N]`` (0 H) and an ammine ``[NH3]`` (3 H) bind the metal as a
lone nitrogen with no heavy neighbour, and before this fix both serialized to the
same bare ``N{n}`` token: ``generate_inline_string``'s ``replace_map`` de-bracketed
the nitride and had a hard-coded ``NH3 -> N{n}`` case. The generator then read bare
``N``, RDKit filled 3 implicit H, and nitrides (WAYHOW/AFIZEV/FUVNER/OPUYAB) came back
as NH3 and failed on atom count.

The fix (encoder-only): emit a bracketed ``[N]{n}`` for the 0-H terminal nitride while
leaving the ammine as bare ``N{n}``. The generator's ``NoImplicit`` guard treats the
bracket atom as authoritative (keeps it 0-H) and the comparator keys ``[N]`` distinct
from an ammine -- both without any decode-side change. The gate is ``GetDegree() == 0``
so a heavy>=1 amido/imido N (S1's exact 0-H convention) stays de-bracketed.
"""

import unittest

from oinsmiles.oin.compare import canonical_roundtrip_key as key
from oinsmiles.oin.inline import OINInlineHandler


def enc(v2):
    return OINInlineHandler.generate_inline_string(v2)


class TestNitrideAmmineEncoding(unittest.TestCase):
    def test_terminal_nitride_is_bracketed(self):
        """A 0-H terminal nitride donor serializes as ``[N]{n}`` (marker), not the
        bare ``N{n}`` an ammine uses."""
        out = enc("[W].[N].[NH3] |g:OCT|w:1.0:0;2.0:1|")
        self.assertIn("[N]{0}", out)
        # the ammine (fragment 2, slot 1) stays bare
        self.assertIn("N{1}", out)
        self.assertNotIn("[N]{1}", out)

    def test_ammine_stays_bare(self):
        """A heavy==0 ammine stays bare ``N{n}`` -- the four passing ammines must be
        untouched by the nitride marker."""
        out = enc("[Cr].[NH3].[NH3] |g:TET|w:1.0:0;2.0:1|")
        self.assertNotIn("[N]{", out)
        self.assertIn("N{0}", out)
        self.assertIn("N{1}", out)

    def test_heavy_amido_stays_bare(self):
        """A 0-H amido/imido N with a heavy neighbour is bare ``N{n}`` by S1's exact
        convention -- the ``GetDegree() == 0`` gate must not bracket it."""
        self.assertIn("N{0}", enc("[W].C[N]C |g:OCT|w:1.1:0|"))  # dimethylamido, heavy==2
        self.assertNotIn("[N]{", enc("[W].C[N]C |g:OCT|w:1.1:0|"))
        self.assertNotIn("[N]{", enc("[Ti].[N]c1ccccc1 |g:OCT|w:1.0:0|"))  # anilido, heavy==1


class TestNitrideAmmineComparator(unittest.TestCase):
    def test_nitride_keys_distinct_from_ammine(self):
        """A bracketed nitride and a bare ammine are different molecules (3 H apart)
        and must NOT collapse to one key."""
        nitride = "[Cr_TET].Cc1cn{0}c2ccc(F)cc12.CC(C)N{1}C(C)C.CC(C)N{2}C(C)C.[N]{3}"
        ammine = "[Cr_TET].Cc1cn{0}c2ccc(F)cc12.CC(C)N{1}C(C)C.CC(C)N{2}C(C)C.N{3}"
        self.assertNotEqual(key(nitride), key(ammine))

    def test_nitride_roundtrips(self):
        """With the marker on both the input encode and the generated re-encode, a
        0-H nitride keys equal to itself (the round trip the fix restores)."""
        nitride = "[W_OCT].CO{0}CCO{3}C.[N]{1}.CC(O{2})(C(F)(F)F)C(F)(F)F"
        self.assertEqual(key(nitride), key(nitride))

    def test_ammine_roundtrips(self):
        """A bare-``N`` ammine still keys equal to itself (no regression)."""
        ammine = "[W_OCT].CO{0}CCO{3}C.N{1}.CC(O{2})(C(F)(F)F)C(F)(F)F"
        self.assertEqual(key(ammine), key(ammine))


if __name__ == "__main__":
    unittest.main()
