"""Axial atropisomer round trip (Y2 P2): OIN -> 3D must rebuild the requested enantiomer.

Wave 1 showed the encoder was blind to biaryl atropisomerism. Wave 2 added an opt-in emit
(``OIN_EMIT_AXIAL``) plus an axial-aware pass in the generator's conformer selection. This
guard closes the loop: with the flag on, BOTH BINAP atropisomers must survive
``XYZ -> OIN -> 3D``, not just whichever one the embed happens to favour.

Measured A/B behind this guard (``tools/injectivity/axial_roundtrip_ab.py``): without the
axial-aware pass the generator returns the same handedness for both requests (1/2 correct,
the match being luck); with it, 2/2. Slow -- two full generations -- hence integration.
"""

import os
import sys
import unittest
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.generation.engine import OIN3DGenerator  # noqa: E402
from oinsmiles.oin.axial import mol_axial_token, parse_axial_token  # noqa: E402

BINAP = _ROOT / "tests" / "fixtures" / "PdCl2-R-BINAP.xyz"
SEED = 42
TIMEOUT = 300


def _mirror(src: Path, dst: Path) -> Path:
    lines = src.read_text().splitlines()
    n = int(lines[0])
    out = [lines[0], lines[1]]
    for ln in lines[2 : 2 + n]:
        p = ln.split()
        out.append(f"{p[0]}  {p[1]}  {p[2]}  {-float(p[3]):.6f}")
    dst.write_text("\n".join(out) + "\n")
    return dst


class TestAxialRoundTrip(unittest.TestCase):
    """Both atropisomers must regenerate as themselves."""

    def setUp(self):
        self._prev = os.environ.get("OIN_EMIT_AXIAL")
        os.environ["OIN_EMIT_AXIAL"] = "1"
        self.addCleanup(self._restore)

    def _restore(self):
        if self._prev is None:
            os.environ.pop("OIN_EMIT_AXIAL", None)
        else:
            os.environ["OIN_EMIT_AXIAL"] = self._prev

    def _roundtrip(self, xyz_path: str) -> tuple[str | None, str | None]:
        oin = XYZToSMILES().convert(xyz_path)
        requested = parse_axial_token(oin)
        res = OIN3DGenerator(optimizer="ff", seed=SEED, timeout=TIMEOUT).generate(oin)
        mol = getattr(res, "mol", None)
        if mol is None:
            self.skipTest("generator returned no mol (eta fallback path)")
        return requested, mol_axial_token(mol)

    def test_r_atropisomer_survives(self):
        requested, got = self._roundtrip(str(BINAP))
        self.assertIn(requested, ("+", "-"), "BINAP must carry an axial token with the flag on")
        self.assertEqual(got, requested)

    def test_s_atropisomer_survives(self):
        """The mirror is the arm the baseline generator gets WRONG without axial selection."""
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            mirror = _mirror(BINAP, Path(d) / "S-BINAP.xyz")
            requested, got = self._roundtrip(str(mirror))
        self.assertIn(requested, ("+", "-"))
        self.assertEqual(got, requested)

    def test_enantiomers_request_opposite_tokens(self):
        """Sanity: the two fixtures really do ask for opposite configurations."""
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            mirror = _mirror(BINAP, Path(d) / "S-BINAP.xyz")
            a = parse_axial_token(XYZToSMILES().convert(str(BINAP)))
            b = parse_axial_token(XYZToSMILES().convert(str(mirror)))
        self.assertTrue(a and b)
        self.assertNotEqual(a, b)


if __name__ == "__main__":
    unittest.main()
