"""Opt-in axial / atropisomer emit (Y2 P2 fix candidate) -- guard the contract.

The encoder is blind to biaryl atropisomerism (Wave 1: R-BINAP and S-BINAP are byte-
identical). ``OIN_EMIT_AXIAL`` opts into an axial-sign token that makes the raw string
distinguish the two enantiomers. The contract these guards lock in:

* **default OFF** -> byte-identical to before (no regression, blindness preserved);
* **flag ON** -> the raw strings of the R and S atropisomers diverge (``|ax:-|`` vs
  ``|ax:+|``), so the encoder is injective on this axis in the raw string;
* the round-trip **key still folds** the token, so the batch harness is unaffected whether
  or not the flag is set (the generator cannot yet reproduce the axis -- emitting a gated
  token would otherwise convert a silent false-positive into a false-negative);
* molecules with **no hindered biaryl axis** are byte-identical even with the flag ON.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key  # noqa: E402

FIX = _ROOT / "tests" / "fixtures"
BINAP = str(FIX / "PdCl2-R-BINAP.xyz")
CISPLATIN = str(FIX / "CisPlatin.xyz")


def _mirror_xyz(path: str, dst: Path) -> str:
    """Write the z-mirror (enantiomer) of an XYZ fixture; return its path."""
    lines = Path(path).read_text().splitlines()
    out = [lines[0], lines[1]]
    for ln in lines[2:]:
        parts = ln.split()
        if len(parts) >= 4:
            parts[3] = f"{-float(parts[3]):.6f}"
            out.append("  ".join(parts))
        elif ln.strip():
            out.append(ln)
    dst.write_text("\n".join(out) + "\n")
    return str(dst)


class _AxialEmitBase(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("OIN_EMIT_AXIAL")
        self.addCleanup(self._restore)
        self._tmp = _ROOT / "tests" / "fixtures"

    def _restore(self):
        if self._prev is None:
            os.environ.pop("OIN_EMIT_AXIAL", None)
        else:
            os.environ["OIN_EMIT_AXIAL"] = self._prev

    def _set(self, on: bool):
        if on:
            os.environ["OIN_EMIT_AXIAL"] = "1"
        else:
            os.environ.pop("OIN_EMIT_AXIAL", None)


class TestDefaultOff(_AxialEmitBase):
    def test_binap_blind_by_default(self):
        # unchanged Wave-1 behaviour: the atropisomers are byte-identical with the flag off.
        self._set(False)

        with tempfile.TemporaryDirectory() as d:
            mir = _mirror_xyz(BINAP, Path(d) / "S.xyz")
            r = XYZToSMILES().convert(BINAP)
            s = XYZToSMILES().convert(mir)
        self.assertEqual(r, s, "default OFF must stay byte-identical (no regression)")
        self.assertNotIn("|ax:", r)


class TestFlagOn(_AxialEmitBase):
    def test_binap_atropisomers_diverge_in_raw_string(self):
        self._set(True)

        with tempfile.TemporaryDirectory() as d:
            mir = _mirror_xyz(BINAP, Path(d) / "S.xyz")
            r = XYZToSMILES().convert(BINAP)
            s = XYZToSMILES().convert(mir)
        self.assertNotEqual(r, s, "flag ON must distinguish R and S in the raw string")
        self.assertTrue(r.endswith("|ax:-|") or r.endswith("|ax:+|"))
        self.assertTrue(s.endswith("|ax:-|") or s.endswith("|ax:+|"))
        self.assertNotEqual(r[-4:], s[-4:], "the two signs must be opposite")

    def test_key_still_folds_the_token(self):
        # the batch harness gates on the key; it must be blind to the opt-in token.
        self._set(True)

        with tempfile.TemporaryDirectory() as d:
            mir = _mirror_xyz(BINAP, Path(d) / "S.xyz")
            kr = winding_canonical_key(normalize_oin_for_comparison(XYZToSMILES().convert(BINAP)))
            ks = winding_canonical_key(normalize_oin_for_comparison(XYZToSMILES().convert(mir)))
        self.assertEqual(kr, ks, "the round-trip key must fold the axial token")

    def test_non_atropisomer_unaffected(self):
        # cisplatin has no hindered biaryl axis: byte-identical even with the flag ON.
        self._set(False)
        off = XYZToSMILES().convert(CISPLATIN)
        self._set(True)
        on = XYZToSMILES().convert(CISPLATIN)
        self.assertEqual(off, on)
        self.assertNotIn("|ax:", on)


class TestTokenIsCanonical(unittest.TestCase):
    """The axial token is a genuine chirality descriptor, not merely a distinguisher.

    A canonical descriptor must depend only on the molecular graph plus the *handedness*
    of the geometry: invariant under input atom renumbering and under any proper rotation,
    flipping only under reflection.
    """

    @staticmethod
    def _atoms(path):

        lines = Path(path).read_text().splitlines()
        n = int(lines[0])
        return lines[1], [
            (ln.split()[0], np.array([float(x) for x in ln.split()[1:4]]))
            for ln in lines[2 : 2 + n]
        ]

    @staticmethod
    def _token(head, atoms):

        from oinsmiles.oin.axial import axial_token
        from tools.injectivity.config_oracle import load_mol

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "probe.xyz"
            body = "".join(f"{s} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}\n" for s, c in atoms)
            p.write_text(f"{len(atoms)}\n{head}\n{body}")
            return axial_token(load_mol(str(p)))

    def test_invariant_under_atom_renumbering(self):

        head, atoms = self._atoms(BINAP)
        base = self._token(head, atoms)
        self.assertTrue(base, "BINAP must produce a token at all")
        rng = np.random.default_rng(42)
        for _ in range(2):
            perm = [0] + [i for i in rng.permutation(len(atoms)) if i != 0]  # metal stays first
            self.assertEqual(self._token(head, [atoms[i] for i in perm]), base)

    def test_invariant_under_proper_rotation(self):

        head, atoms = self._atoms(BINAP)
        base = self._token(head, atoms)
        rng = np.random.default_rng(7)
        for _ in range(2):
            q, r = np.linalg.qr(rng.normal(size=(3, 3)))
            q = q @ np.diag(np.sign(np.diag(r)))
            if np.linalg.det(q) < 0:
                q[:, 0] *= -1  # force a PROPER rotation
            self.assertEqual(self._token(head, [(s, q @ c) for s, c in atoms]), base)

    def test_flips_under_reflection(self):

        head, atoms = self._atoms(BINAP)
        base = self._token(head, atoms)
        mirrored = self._token(head, [(s, np.array([c[0], c[1], -c[2]])) for s, c in atoms])
        self.assertTrue(mirrored)
        self.assertNotEqual(mirrored, base)


class TestSymmetryEquivalentAxesStillFlip(unittest.TestCase):
    """Regression: a multi-axis token must not be sorted BY SIGN.

    Symmetry-equivalent axes tie on symmetry rank. An earlier version broke that tie with
    the sign itself, to keep the string order-independent -- which silently forced the signs
    into ascending order, so a molecule carrying ``+-`` rendered identically to its mirror
    carrying ``-+``. The token stopped being a chirality descriptor for exactly the
    structures that need it most.

    `YESKOZ` is a real corpus example: two symmetry-equivalent hindered biaryl axes of
    opposite sign, which the independent geometric oracle calls chiral (mirror RMSD 3.2 A).
    Found by the corpus-wide sign-convention audit (34/37 flipped; this was one of the 3).
    """

    def test_two_equivalent_axes_of_opposite_sign_flip(self):

        from oinsmiles.oin.axial import axial_token
        from tools.injectivity.config_oracle import load_mol

        fixture = FIX / "YESKOZ.xyz"
        base = axial_token(load_mol(str(fixture)))
        self.assertEqual(len(base), 2, "YESKOZ must present exactly two emitting axes")
        self.assertEqual(set(base), {"+", "-"}, "the two signs must be opposite")

        lines = fixture.read_text().splitlines()
        n = int(lines[0])

        with tempfile.TemporaryDirectory() as d:
            mirrored = Path(d) / "mirror.xyz"
            body = []
            for ln in lines[2 : 2 + n]:
                p = ln.split()
                body.append(f"{p[0]}  {p[1]}  {p[2]}  {-float(p[3]):.6f}")
            mirrored.write_text(f"{lines[0]}\n{lines[1]}\n" + "\n".join(body) + "\n")
            mirror = axial_token(load_mol(str(mirrored)))

        expected = base.translate(str.maketrans("+-", "-+"))
        self.assertEqual(mirror, expected, "the multi-axis token must flip for the mirror")
        self.assertNotEqual(mirror, base)


class TestNotOverSensitive(unittest.TestCase):
    """A hindered biaryl with a symmetry-degenerate ring end is ACHIRAL -- no token.

    A ring end whose two ortho neighbours are symmetry-equivalent carries a local mirror
    plane through the axis, so the molecule is not an atropisomer however twisted it is.
    Emitting a sign there would make the encoder *over-sensitive*: claiming a stereo
    distinction that does not exist (the `over_sensitive` cell of the Y1 confusion matrix).
    """

    @staticmethod
    def _from_smiles(smi):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)
        return m

    def test_symmetric_ring_end_emits_nothing(self):
        from oinsmiles.oin.axial import axial_token, detect_axial_axes

        # 2,6-dimethyl on one end -> identical ortho substituents -> achiral
        m = self._from_smiles("Cc1cccc(C)c1-c1c(C)cccc1OC")
        hindered = [a for a in detect_axial_axes(m) if a.hindered]
        self.assertTrue(hindered, "the axis IS sterically hindered ...")
        self.assertFalse(any(a.stereogenic for a in hindered), "... but it is NOT stereogenic")
        self.assertEqual(axial_token(m), "", "an achiral biaryl must emit no axial token")

    def test_asymmetric_ring_ends_do_emit(self):
        from oinsmiles.oin.axial import axial_token

        m = self._from_smiles("COc1cccc(C)c1-c1c(C)cccc1[N+](=O)[O-]")
        self.assertIn(axial_token(m), ("+", "-"))

    def test_gate_agrees_with_independent_oracle(self):
        """Cross-validate the gate against the geometric oracle, which never sees its logic.

        These synthetic biaryls are rigid, so the oracle (min proper-rotation mirror RMSD) is
        valid for them -- unlike flexible dataset conformers, where it over-reports.
        """
        from oinsmiles.oin.axial import axial_token
        from tools.injectivity.oracle import geometric_chirality

        for smi, expect_chiral in (
            ("COc1cccc(C)c1-c1c(C)cccc1[N+](=O)[O-]", True),  # both ends asymmetric
            ("Cc1cccc(C)c1-c1c(C)cccc1OC", False),  # one end 2,6-dimethyl -> achiral
            ("Cc1cccc(C)c1-c1c(C)cccc1C", False),  # both ends 2,6-dimethyl -> achiral
        ):
            with self.subTest(smiles=smi):
                m = self._from_smiles(smi)
                _, _, chiral = geometric_chirality(m, m.GetConformer().GetPositions())
                self.assertEqual(chiral, expect_chiral, "oracle disagrees with the premise")
                self.assertEqual(
                    bool(axial_token(m)), chiral, "the gate must agree with the oracle"
                )


class TestNormalizerFolds(_AxialEmitBase):
    def test_normalize_strips_axial_suffix(self):
        base = "[Pd_SPL].c1ccccc1{0}.[Cl]{1}"
        self.assertEqual(
            normalize_oin_for_comparison(base + " |ax:-|"),
            normalize_oin_for_comparison(base),
        )


if __name__ == "__main__":
    unittest.main()
