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

    def test_emit_gate_is_off_unless_the_env_var_is_set(self):
        """``OIN_EMIT_AXIAL`` must stay **opt-in**. Guard against an accidental promotion.

        v0.4.5 product call: injectivity levers stay default-OFF for this release. The
        evidence to promote exists (cohort A/B 22/22 vs 8/22 on single-axis structures) and is
        written up as a staged v0.4.6 recommendation in ``docs/KNOWN_LIMITATIONS.md``, but the
        flip is deliberately not taken here -- turning it on changes the emitted string for
        part of the corpus and makes the round-trip key's fold of the token load-bearing.

        Unlike ``OIN_EARLY_EXIT`` (``os.environ.get(..., "1") != "0"``) this gate is a bare
        truthiness test on an unset variable, so "off" is the absence of the variable. Assert
        that shape directly: an unset environment must not emit.
        """
        self._set(False)
        self.assertIsNone(os.environ.get("OIN_EMIT_AXIAL"))

        oin = XYZToSMILES().convert(BINAP)
        self.assertNotIn("|ax:", oin, "OIN_EMIT_AXIAL must default to OFF in v0.4.5")

        os.environ["OIN_EMIT_AXIAL"] = "0"
        self.assertNotIn(
            "|ax:", XYZToSMILES().convert(BINAP), "an explicit falsey value must stay OFF"
        )


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


def _two_axis_probe(twists=(75.0, -75.0)):
    """A molecule with two SYMMETRY-EQUIVALENT hindered axes, twisted as asked.

    Two disconnected copies of 2,2'-dimethylbiphenyl. Each axis end is a 2-methylphenyl,
    whose ortho pair (C-CH3 vs C-H) is asymmetric, so both axes are hindered and stereogenic;
    the two copies are indistinguishable, so the axes *tie* on symmetry rank -- which is
    precisely the configuration the sign-sorting bug corrupted. Built here rather than taken
    from a fixture so the two signs can be set independently and the guard cannot be
    satisfied by a molecule whose axes happen to agree.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    from oinsmiles.oin.axial import detect_axial_axes

    mol = Chem.AddHs(Chem.MolFromSmiles("Cc1ccccc1-c1ccccc1C.Cc1ccccc1-c1ccccc1C"))
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(mol, ps) != 0:
        raise unittest.SkipTest("ETKDG could not embed the two-axis probe")
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    for ax, want in zip(detect_axial_axes(mol), twists):
        n1 = next(
            n.GetIdx()
            for n in mol.GetAtomWithIdx(ax.a1).GetNeighbors()
            if n.GetIdx() != ax.a2 and n.GetIsAromatic()
        )
        n2 = next(
            n.GetIdx()
            for n in mol.GetAtomWithIdx(ax.a2).GetNeighbors()
            if n.GetIdx() != ax.a1 and n.GetIsAromatic()
        )
        rdMolTransforms.SetDihedralDeg(conf, n1, ax.a1, ax.a2, n2, want)
    return mol


def _z_mirror(mol):
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, -p.z))
    return out


class TestSymmetryEquivalentAxesStillFlip(unittest.TestCase):
    """Regression: a multi-axis token must not be sorted BY SIGN.

    Symmetry-equivalent axes tie on symmetry rank. An earlier version broke that tie with
    the sign itself, to keep the string order-independent -- which silently forced the signs
    into ascending order, so a molecule carrying ``+-`` rendered identically to its mirror
    carrying ``-+``. The token stopped being a chirality descriptor for exactly the
    structures that need it most.

    This used to be guarded with the ``YESKOZ`` fixture, which no longer emits: its two
    meso-aryl axes are **not** stereogenic per-axis (see
    ``TestPorphyrinMesoAxesAreNotPerAxisStereogenic``). A constructed probe is a better guard
    anyway -- the two signs can be set independently, so a token that collapsed them could
    not pass by luck.
    """

    def test_two_equivalent_axes_of_opposite_sign_flip(self):
        from oinsmiles.oin.axial import axial_token

        mol = _two_axis_probe()
        base = axial_token(mol)
        self.assertEqual(len(base), 2, "the probe must present exactly two emitting axes")
        self.assertEqual(set(base), {"+", "-"}, "the two signs must be opposite")

        mirror = axial_token(_z_mirror(mol))
        self.assertEqual(
            mirror,
            base.translate(str.maketrans("+-", "-+")),
            "the multi-axis token must flip for the mirror",
        )
        self.assertNotEqual(mirror, base)

    def test_same_sign_pair_is_not_confused_with_opposite_sign_pair(self):
        """``++`` and ``+-`` must be different strings -- the collapse the old sort caused."""
        from oinsmiles.oin.axial import axial_token

        same = axial_token(_two_axis_probe(twists=(75.0, 75.0)))
        opposite = axial_token(_two_axis_probe(twists=(75.0, -75.0)))
        self.assertEqual(len(same), 2)
        self.assertEqual(len(set(same)), 1, "both axes were twisted the same way")
        self.assertNotEqual(same, opposite)


class TestPerceptionInvariance(unittest.TestCase):
    """The token must not depend on which bond-order model perceived the molecule.

    The encoder perceives bond orders from interatomic distances (``xyz2mol``); the generator
    transfers them from the OIN fragment SMILES (``build_contract_mol``). The two disagree
    for a metalloporphyrin -- the encoder reads an aromatic pyrrolide core on Zn(II), the
    generated mol a neutral localized tautomer. A descriptor keyed on ``GetIsAromatic()``
    therefore found a meso-aryl axis on the input and **none at all** on any generated
    conformer, which is what made the multi-axis cohort unreproducible: the generator holds
    both hindered twists, but nothing downstream could see them.

    The perturbation below changes only the perception -- coordinates, elements and
    connectivity are untouched, so the handedness is untouched -- and the token must be
    unchanged. ``tools/injectivity/axial_perception_sweep.py`` runs the same check at corpus
    scale.
    """

    def _delocalized(self, mol):
        from tools.injectivity.axial_perception_sweep import delocalize

        return delocalize(mol)

    def test_binap_token_survives_delocalization(self):
        from oinsmiles.oin.axial import axial_token
        from tools.injectivity.config_oracle import load_mol

        mol = load_mol(BINAP)
        base = axial_token(mol)
        self.assertTrue(base, "BINAP must emit a token to make this guard meaningful")
        self.assertEqual(axial_token(self._delocalized(mol)), base)

    def test_two_axis_probe_survives_delocalization(self):
        from oinsmiles.oin.axial import axial_token

        mol = _two_axis_probe()
        base = axial_token(mol)
        self.assertEqual(len(base), 2)
        self.assertEqual(axial_token(self._delocalized(mol)), base)


class TestPorphyrinMesoAxesAreNotPerAxisStereogenic(unittest.TestCase):
    """``YESKOZ`` must NOT emit -- its meso-aryl axes have a local C2 at the porphyrin end.

    A 5,15-diarylporphyrin's meso carbon is flanked by two pyrrole rings that a graph
    automorphism swaps, so rotating that half 180 deg about the meso-aryl axis reproduces the
    molecule and the axis carries no per-axis handedness. Both configurations are in fact
    achiral: for the anti (alpha/beta) isomer the mirror plane through the two meso positions
    composed with the C2 about the porphyrin normal maps it to itself, and for syn the mirror
    plane through the other two meso positions does.

    The old descriptor called these axes stereogenic, but only by accident: it ranked atoms
    on the molecule *as perceived*, and the arbitrary resonance form ``AC2BO`` happened to
    return broke the macrocycle's symmetry. So the emitted ``+-`` was a false positive
    resting on a non-canonical choice -- which is why the generator could never reproduce it,
    and why the multi-axis cohort measured 0/2 in *both* A/B arms.

    What is genuinely lost is the syn/anti *diastereomerism* (which face each ortho
    substituent points to). That is a relative configuration across two axes, not a sign on
    one, and a per-axis signed dihedral cannot express it -- see ``docs/KNOWN_LIMITATIONS.md``.
    """

    def test_yeskoz_emits_no_token(self):
        from oinsmiles.oin.axial import axial_token, detect_axial_axes
        from tools.injectivity.config_oracle import load_mol

        mol = load_mol(str(FIX / "YESKOZ.xyz"))
        axes = detect_axial_axes(mol)
        self.assertEqual(len(axes), 2, "the two meso-aryl axes must still be DETECTED")
        self.assertTrue(all(a.hindered for a in axes), "both are sterically hindered")
        self.assertFalse(
            any(a.stereogenic for a in axes),
            "neither meso-aryl axis is stereogenic on its own",
        )
        self.assertEqual(axial_token(mol), "")


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
