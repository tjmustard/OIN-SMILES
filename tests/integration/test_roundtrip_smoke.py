"""Round-trip accuracy smoke test (CI guard).

Regenerates 3D structures for the canonical golden OINs and re-encodes them,
asserting the OIN round-trips (OIN -> XYZ -> OIN) under the deterministic
``optimizer="ff"`` + ``seed=42`` path. This pins the headline round-trip
accuracy *in CI*, so a generation- or encoding-side regression fails the build
instead of only surfacing in the out-of-band dataset harness
(``tools/test_dataset_roundtrip.py``, which is not run in CI).

Design:
- **Deterministic** -- fixed seed + FF optimizer (no g-xTB/MACE binaries), so it
  is hermetic and reproducible.
- **Comparison** reuses the same canonical comparator the full round-trip
  harness uses (``winding_canonical_key(normalize_oin_for_comparison(...))``),
  and the same re-encode path (``get_oin_string`` on the generator's bonded mol).
- **Fixtures** -- the four perf-wave A/B goldens (square-planar, haptic/winding,
  octahedral aromatic chelate, chelating P-donors + axial), plus transplatin so
  a cis/trans slot-ordering flip is caught (cisplatin alone would not).

Runtime ~50s per interpreter (dominated by fac-Ir(ppy)3 and PdCl2-BINAP).
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key
from oinsmiles.utils.perception_tmc import get_oin_string

SEED = 42

# Canonical golden OINs. The first four are the perf-wave A/B set; transplatin
# guards cis/trans slot-ordering (cisplatin alone would not catch an isomer flip).
GOLDEN_OINS = {
    "cisplatin": "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
    "transplatin": "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}",
    "ferrocene": (
        "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    ),
    "fac-Ir(ppy)3": (
        "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
    ),
    "PdCl2-BINAP": (
        "[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}(c3ccccc3)c3ccccc3)"
        "ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}"
    ),
}


def _roundtrip_key(oin_string):
    """Canonical isomer key used for OIN string identity (matches the harness)."""
    return winding_canonical_key(normalize_oin_for_comparison(oin_string.strip()))


def _reencode(gen_mol, xyz):
    """Re-encode a generated structure to OIN via the generator's bonded mol.

    Mirrors tests/integration/verify_roundtrip.py: parse coords from the XYZ
    block and hand the bonded mol to get_oin_string.
    """
    lines = xyz.splitlines()
    natoms = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + natoms)])
    return get_oin_string(gen_mol, coords)


class TestRoundTripSmoke(unittest.TestCase):
    """OIN -> XYZ -> OIN string identity on the golden fixtures (FF + fixed seed)."""

    def test_goldens_roundtrip(self):
        for name, oin1 in GOLDEN_OINS.items():
            with self.subTest(fixture=name):
                gen = OIN3DGenerator(optimizer="ff", seed=SEED).generate(oin1)
                self.assertIsNotNone(
                    gen.mol,
                    f"{name}: generator returned no bonded mol; cannot re-encode",
                )
                oin2 = _reencode(gen.mol, gen.xyz)
                self.assertEqual(
                    _roundtrip_key(oin1),
                    _roundtrip_key(oin2),
                    f"{name}: OIN round-trip mismatch\n"
                    f"  in : {normalize_oin_for_comparison(oin1)}\n"
                    f"  out: {normalize_oin_for_comparison(oin2)}",
                )


if __name__ == "__main__":
    unittest.main()
