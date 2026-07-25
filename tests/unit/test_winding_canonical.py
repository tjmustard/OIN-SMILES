"""Lane 3 guards: the eta winding marker is orientation-invariant, and canonicalizing
it never folds an enantiomer.

Two independent properties are pinned here.

**Orientation invariance** (the lane's definition of done). Re-encoding a structure from
a random *proper* rotation must give a byte-identical OIN string. That covers both halves
of the marker at once -- a moved heading atom and a flipped ``>``/``<`` both change the
string -- so it is a stronger assertion than checking either alone. No XYZ hash is
hardcoded: a prior wave flaked CI that way, so every assertion here compares the encoder
against *itself* under a transform.

**Chirality preservation** (the anti-trap guard). `OIN_CANONICAL_ETA_WINDING` collapses
the two mirror-related spellings of an *achiral* eta arrangement. The v0.4.4 axial wave
shipped a canonicalization that looked right on one easy fixture and was in fact
reflection-invariant -- it destroyed the stereochemistry it encoded. So the guard here is
not "the fixture still encodes": it is that `TiCat3` and `TiCat4`, a real
**rac / meso pair of the same Me2Si-bridged bis(indenyl) ligand**, stay distinct with the
lever on, and that a chiral structure still encodes differently from its mirror image.

`TiCat3` / `TiCat4` are load-bearing precisely because they differ *only* in eta winding:
same metal, same geometry, same ligand bodies. If a canonicalization ever folds them, the
notation has stopped recording which face of each ring the metal is on.
"""

import os
import unittest
from pathlib import Path

import numpy as np

from oinsmiles import XYZToSMILES
from oinsmiles.utils import oin_aligner as AL

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"

#: Eta fixtures spanning the awkward cases, not just the easy one: an unsubstituted
#: metallocene (orientation-free rings, winding forced), a substituted-face ferrocene,
#: a bridged bis(Cp) constrained-geometry catalyst, and the rac/meso bis(indenyl) pair.
ETA_FIXTURES = [
    "Ferrocene.xyz",
    "Ferrocene-halide-face.xyz",
    "TiCat1.xyz",
    "TiCat3.xyz",
    "TiCat4.xyz",
]


def _read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for line in lines[2 : 2 + n]:
        p = line.split()
        syms.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.array(xyz), (lines[1] if len(lines) > 1 else "")


def _write_xyz(path, syms, xyz, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, xyz):
            f.write(f"{s} {x:.10f} {y:.10f} {z:.10f}\n")


def _random_proper_rotation(rng):
    """Uniform random rotation with det = +1 -- never a reflection."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    assert np.linalg.det(q) > 0
    return q


class _EtaEncodeCase(unittest.TestCase):
    def setUp(self):
        self.conv = XYZToSMILES()
        self.tmp = Path(os.environ.get("TMPDIR", "/tmp")) / "oin_winding_canonical_test"
        self.tmp.mkdir(parents=True, exist_ok=True)

    def encode(self, path):
        return self.conv.convert(str(path))

    def encode_transformed(self, path, matrix, tag):
        syms, xyz, comment = _read_xyz(path)
        out = self.tmp / f"{Path(path).stem}_{tag}.xyz"
        _write_xyz(out, syms, xyz @ matrix.T, comment)
        return self.encode(out)


class TestWindingOrientationInvariance(_EtaEncodeCase):
    """The heading atom AND the sign survive a random rigid rotation of the input."""

    def test_eta_encoding_is_invariant_under_random_proper_rotations(self):
        rng = np.random.default_rng(20450103)
        for name in ETA_FIXTURES:
            path = FIXTURES / name
            with self.subTest(fixture=name):
                base = self.encode(path)
                self.assertIsNotNone(base)
                for trial in range(3):
                    rotated = self.encode_transformed(
                        path, _random_proper_rotation(rng), f"rot{trial}"
                    )
                    self.assertEqual(
                        base,
                        rotated,
                        f"{name}: OIN changed under a proper rotation "
                        f"(heading atom or winding sign is embedding-dependent)",
                    )


class TestEtaWindingCanonicalizationPreservesChirality(_EtaEncodeCase):
    """The gated canonicalization must never make the encoder reflection-invariant."""

    @staticmethod
    def _mirror(xyz):
        out = xyz.copy()
        out[:, 0] *= -1.0
        return out

    def _encode_with_lever(self, path, enabled):
        original = AL.CANONICAL_ETA_WINDING
        AL.CANONICAL_ETA_WINDING = enabled
        try:
            return self.encode(path)
        finally:
            AL.CANONICAL_ETA_WINDING = original

    def _encode_mirror_with_lever(self, path, enabled):
        syms, xyz, comment = _read_xyz(path)
        out = self.tmp / f"{Path(path).stem}_mirror.xyz"
        _write_xyz(out, syms, self._mirror(xyz), comment)
        return self._encode_with_lever(out, enabled)

    def test_rac_and_meso_bis_indenyl_stay_distinct_with_lever_on(self):
        # TiCat3 and TiCat4 are the SAME ligand on the SAME metal in the SAME geometry;
        # only the coordinated ring faces differ. Folding them would erase rac/meso.
        for enabled in (False, True):
            with self.subTest(lever=enabled):
                three = self._encode_with_lever(FIXTURES / "TiCat3.xyz", enabled)
                four = self._encode_with_lever(FIXTURES / "TiCat4.xyz", enabled)
                self.assertNotEqual(
                    three,
                    four,
                    "rac and meso bis(indenyl) collapsed to one string -- the winding "
                    "marker has stopped recording the coordinated face",
                )

    def test_chiral_eta_structure_still_differs_from_its_mirror_with_lever_on(self):
        # TiCat3 is the chiral member of the pair (confirmed independently by the
        # rigid-core mirror test in tools/eta_core_chirality.py: CHIRAL, 2.14 A, vs
        # TiCat4 ACHIRAL, 0.03 A). Its mirror is a genuinely different isomer and must
        # keep a different string whether or not the lever is on.
        for enabled in (False, True):
            with self.subTest(lever=enabled):
                base = self._encode_with_lever(FIXTURES / "TiCat3.xyz", enabled)
                mirror = self._encode_mirror_with_lever(FIXTURES / "TiCat3.xyz", enabled)
                self.assertNotEqual(
                    base,
                    mirror,
                    "a chiral eta complex encodes identically to its mirror image -- "
                    "the canonicalization is reflection-invariant (the v0.4.4 axial bug)",
                )

    def test_achiral_eta_structure_folds_onto_its_mirror_only_with_lever_on(self):
        # TiCat4 is the achiral member. Its mirror is the SAME compound, so the two
        # spellings must collapse once the lever is on -- and, as the honest cost of a
        # default-OFF gate, must still differ with it off.
        path = FIXTURES / "TiCat4.xyz"
        self.assertNotEqual(
            self._encode_with_lever(path, False),
            self._encode_mirror_with_lever(path, False),
            "precondition: with the lever OFF the achiral pair is expected to drift",
        )
        self.assertEqual(
            self._encode_with_lever(path, True),
            self._encode_mirror_with_lever(path, True),
            "with the lever ON the two mirror-related spellings of an ACHIRAL eta "
            "arrangement must canonicalize to one string",
        )

    def test_lever_leaves_orientation_free_metallocenes_untouched(self):
        # Cp rings can be turned over by a proper rotation, so their winding is notation,
        # already pinned to a fixed character. The lever must not touch them.
        for name in ("Ferrocene.xyz", "TiCat1.xyz"):
            with self.subTest(fixture=name):
                self.assertEqual(
                    self._encode_with_lever(FIXTURES / name, False),
                    self._encode_with_lever(FIXTURES / name, True),
                )


class TestEtaSwapSense(unittest.TestCase):
    """The cyclic-sense factor must be COMPUTED -- it is not a constant."""

    def test_sense_is_measured_per_fragment_not_assumed(self):
        # Getting this factor wrong inverts the achirality test and silently folds
        # enantiomers, which is exactly what a first cut of this code did. Two separate
        # copies of one fragment are canonically ordered identically (+1); two rings
        # inside ONE bridged fragment need the real automorphism, and the answer differs
        # between real ligands.
        separate = AL.OINDiscreteAligner._eta_swap_sense(
            {"smiles": "c1cc[cH]c1", "rank": 1, "cons": [0, 1, 2, 3, 4]},
            {"smiles": "c1cc[cH]c1", "rank": 2, "cons": [0, 1, 2, 3, 4]},
        )
        self.assertEqual(separate, 1)

        # Different ligand bodies are not automorphic at all -> no fold is permitted.
        self.assertIsNone(
            AL.OINDiscreteAligner._eta_swap_sense(
                {"smiles": "c1cc[cH]c1", "rank": 1, "cons": [0, 1, 2, 3, 4]},
                {"smiles": "Cc1cc[cH]c1", "rank": 2, "cons": [1, 2, 3, 4, 5]},
            )
        )


if __name__ == "__main__":
    unittest.main()
