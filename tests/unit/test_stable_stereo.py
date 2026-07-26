"""Correctness and atom-numbering stability of tetrahedral stereo tags (v0.4.5 Lane 8).

The defect
==========
``get_oin_string`` rebuilds each ligand fragment atom-by-atom. ``RWMol.AddAtom``
copies the parent's chiral tag **verbatim**, but a chiral tag is a parity *relative
to the neighbour order on that atom* -- and the rebuild changes that order twice,
both times as a function of the input atom numbering: hydrogens are folded into
``SetNumExplicitHs``, and bonds are re-added by ascending parent index. So the
same 3D structure, presented with its atoms in a different order, can emit ``@``
where it previously emitted ``@@`` (measured at 13.0% of a 225-molecule sample) --
and on any single ordering the emitted configuration is right only by luck.

``OIN_STABLE_STEREO`` (default OFF) re-derives the fragment tags from the parent
geometry instead of translating them.

Two properties, tested together
===============================
Stability alone is worthless: a descriptor that is *constant* is perfectly stable
and encodes nothing. (A Y2-wave canonicalization shipped exactly that failure --
it made an axial descriptor reflection-invariant, and every guard written against
the one easy fixture passed.) So the stability assertions here are paired with:

* :meth:`TestStableStereoIsCorrect.test_mirror_image_inverts_every_tag` -- the
  descriptor must still flip for the enantiomer;
* :meth:`TestStableStereoIsCorrect.test_lever_on_emits_the_true_configuration` --
  it must agree with the absolute configuration the fixture is *named* for; and
* :meth:`TestStableStereoIsCorrect.test_tag_count_is_preserved` -- stability must
  not be bought by emitting fewer tags.
"""

import os
import random
import re
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdCIPLabeler  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))
_LEVER = "OIN_STABLE_STEREO"

# One stereo tag, counted once whether written @ or @@.
_TAG_RE = re.compile(r"@{1,2}")
_SLOT_RE = re.compile(r"\{\d+[<>^]?\}")
_GEO_RE = re.compile(r"_[A-Z]{3}(?=\])")

# 49 atoms, one Rh, tetrahedral carbon stereocentres on a bidentate ligand. Taken
# from the 29 renumbering-unstable molecules found by tools/canonicality_probe.py
# because it is small, encodes quickly, and drifts in the stereo tags WITHOUT any
# accompanying skeleton or slot drift -- so a failure here can only be this defect.
_UNSTABLE = "ROGYAO_comp_0.xyz"

# (2R,4R)-pentane-2,4-diyl backbones. Their absolute configuration is checkable by
# hand from the raw coordinates -- priority at each centre is D > CH2 > CH3 > H
# with D = P or N, and the triple product of the unit vectors to priorities 1,2,3
# is negative for both centres in both fixtures. So the correct answer is (R,R),
# independently of anything RDKit or this codebase computes.
_RR_FIXTURES = ("PdCl2-RR-BDPP.xyz", "PdCl2-RR-BDNN.xyz")


def _fixture(name):
    return os.path.join(_FIXTURES, name)


def _read_xyz(path):
    with open(path) as fh:
        lines = fh.readlines()
    n = int(lines[0].split()[0])
    syms, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, coords


def _write_xyz(path, syms, coords):
    with open(path, "w") as fh:
        fh.write(f"{len(syms)}\n\n")
        for sym, (x, y, z) in zip(syms, coords):
            fh.write(f"{sym:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")


def _renumbered(path, tmpdir, seed):
    """The same structure with its atoms written in a different order.

    Coordinates are carried across untouched, so the geometry -- and therefore the
    correct encoding -- is bit-for-bit the same molecule.
    """
    syms, coords = _read_xyz(path)
    order = list(range(len(syms)))
    random.Random(seed).shuffle(order)
    out = os.path.join(tmpdir, f"renum_{seed}.xyz")
    _write_xyz(out, [syms[i] for i in order], [coords[i] for i in order])
    return out


def _mirrored(path, tmpdir):
    """Reflected through the yz plane: an improper operation, so the enantiomer."""
    syms, coords = _read_xyz(path)
    out = os.path.join(tmpdir, "mirror.xyz")
    _write_xyz(out, syms, [[-x, y, z] for x, y, z in coords])
    return out


def _encode(path, stable):
    # Set the lever EXPLICITLY in both directions. Unsetting it used to mean "off", but it
    # was promoted to default-ON in v0.4.5 -- at which point `stable=False` silently became
    # `stable=True`, and the two "the defect still reproduces" guards below began asserting
    # the defect against the *fixed* path. They failed loudly, which is the only reason this
    # was caught; a guard that goes vacuous usually just goes green.
    env = dict(os.environ)
    env[_LEVER] = "1" if stable else "0"
    with mock.patch.dict(os.environ, env, clear=True):
        return XYZToSMILES().convert(path)


def _tags(oin):
    return _TAG_RE.findall(oin)


def _skeleton(oin):
    return _TAG_RE.sub("", oin)


def _cip_labels(oin):
    """CIP descriptors of every stereocentre in the emitted string.

    Reads the OIN back the way a consumer would: strip the OIN-specific decorations
    ({slot} markers, the _GEO suffix on the metal) and let ``rdCIPLabeler`` -- the
    rigorous implementation, not the legacy one -- assign R/S.
    """
    labels = []
    for frag in oin.split("|")[0].split("."):
        body = _GEO_RE.sub("", _SLOT_RE.sub("", frag))
        mol = Chem.MolFromSmiles(body)
        if mol is None:
            continue
        if all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
            continue
        rdCIPLabeler.AssignCIPLabels(mol)
        labels += [a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")]
    return labels


class TestStereoIsStableUnderRenumbering(unittest.TestCase):
    """With the lever ON, atom numbering must not change the configuration."""

    # A WIDE seed set, deliberately. The two assertions below are asymmetric: stability
    # must hold for EVERY seed, while "the fixture really is unstable" needs only ONE seed
    # to drift. So widening strengthens the first and de-flakes the second at once.
    #
    # It was originally (1, 2, 3, 5, 8) and that was too narrow: none of those five happens
    # to trigger this molecule's drift, so `test_lever_off_reproduces_the_defect` failed even
    # though the molecule genuinely IS unstable -- `canonicality_probe.py --only
    # ROGYAO_comp_0 --trials 5`, which draws its own permutations, sees 5 of 15
    # transform-trials drift (subclass `rdkit_canonical`, key-level too). A five-sample miss
    # on a ~1-in-3 event is unremarkable; hardcoding five seeds was the bug, not the fixture.
    SEEDS = tuple(range(12))

    def test_lever_off_reproduces_the_defect(self):
        """Guard against a vacuous suite: the fixture must really be unstable.

        If this stops failing, the fixture no longer exercises the bug and the stability
        assertion below has quietly become worthless. Keep this test -- it is what stops
        the lane from proving nothing -- but note it asserts over the whole seed set, not
        any particular seed, because which permutations trigger the drift is not a
        property worth pinning.
        """
        path = _fixture(_UNSTABLE)
        base = _encode(path, stable=False)
        with tempfile.TemporaryDirectory() as tmp:
            drifted = [
                seed
                for seed in self.SEEDS
                if _encode(_renumbered(path, tmp, seed), stable=False) != base
            ]
        self.assertTrue(
            drifted,
            f"{_UNSTABLE} did not drift under any of {len(self.SEEDS)} renumberings with "
            f"{_LEVER} unset. Either the defect was fixed elsewhere (good -- then retire this "
            "fixture and say so) or the fixture no longer exercises it, in which case the "
            "stability assertion below proves nothing. Cross-check with: "
            "tools/canonicality_probe.py --only ROGYAO_comp_0 --trials 5",
        )

    def test_lever_on_is_byte_stable_under_renumbering(self):
        path = _fixture(_UNSTABLE)
        base = _encode(path, stable=True)
        self.assertIn("@", base, "fixture must carry stereochemistry to be a valid test")
        with tempfile.TemporaryDirectory() as tmp:
            for seed in self.SEEDS:
                with self.subTest(seed=seed):
                    self.assertEqual(_encode(_renumbered(path, tmp, seed), stable=True), base)

    def test_lever_on_keeps_reviewed_complexes_stable(self):
        for name in _RR_FIXTURES + ("Rh-RR-DIPAMP-Cl2.xyz",):
            path = _fixture(name)
            base = _encode(path, stable=True)
            with tempfile.TemporaryDirectory() as tmp:
                for seed in (4, 9):
                    with self.subTest(fixture=name, seed=seed):
                        self.assertEqual(_encode(_renumbered(path, tmp, seed), stable=True), base)


class TestStableStereoIsCorrect(unittest.TestCase):
    """Stable is not enough -- the tag has to be the right one, and still chiral."""

    def test_lever_on_emits_the_true_configuration(self):
        """(2R,4R) fixtures must encode as R,R.

        The reference is not another code path: priority at each backbone centre is
        D > CH2 > CH3 > H by inspection, and the triple product of the unit vectors
        to priorities 1,2,3 is negative in the raw coordinates for both centres of
        both fixtures, which is R.
        """
        for name in _RR_FIXTURES:
            with self.subTest(fixture=name):
                self.assertEqual(_cip_labels(_encode(_fixture(name), stable=True)), ["R", "R"])

    def test_lever_off_emits_the_inverted_configuration(self):
        """Documents the defect being fixed: today's encoder reports the enantiomer.

        These two fixtures encode as S,S on the unpatched path -- the wrong absolute
        configuration, not merely an unstable one. If this ever starts returning
        R,R the defect is gone from the default path and this whole lever can be
        promoted and this test deleted.
        """
        for name in _RR_FIXTURES:
            with self.subTest(fixture=name):
                self.assertEqual(_cip_labels(_encode(_fixture(name), stable=False)), ["S", "S"])

    def test_mirror_image_inverts_every_tag(self):
        """Reflect the input: the constitution must survive and every tag must flip.

        A descriptor made stable by being constant passes every stability test above
        and fails this one.
        """
        # _UNSTABLE (ROGYAO_comp_0) is deliberately NOT in this list. It carries two @@ tags
        # and looks like an obvious enantiomer test, but it is ACHIRAL: the corrected
        # torsion-aware oracle superimposes its mirror at 0.423 A with only 4 automorphisms
        # (RIGID_ACHIRAL, threshold 0.5). Demanding that its mirror differ asserts something
        # false. It stays in the renumbering-stability tests above, where its instability is
        # real and measured -- an unstable fixture and a chiral fixture are different jobs.
        #
        # Worth knowing WHY it looked chiral: the older rigid oracle
        # (tools/injectivity/oracle.py) reported "distinct, mirror RMSD 2.586 A,
        # ENCODER-BLIND (total)" -- a false positive from hitting its 4000-automorphism cap.
        # It enumerates on the H-EXPLICIT graph, where this molecule's four methyls starve
        # the budget, which inflates rigid RMSD on methyl-rich species. That instrument
        # defect is recorded in the Lane 7 write-up; it changes no curated fixture verdict
        # but it will send you chasing bugs that do not exist. Use torsion_oracle.py.
        for name in _RR_FIXTURES:
            with self.subTest(fixture=name):
                path = _fixture(name)
                base = _encode(path, stable=True)
                self.assertIn("@", base)
                with tempfile.TemporaryDirectory() as tmp:
                    mirror = _encode(_mirrored(path, tmp), stable=True)
                self.assertEqual(
                    _skeleton(base),
                    _skeleton(mirror),
                    "reflection must not change the constitution",
                )
                self.assertNotEqual(base, mirror, "the enantiomer must encode differently")
                self.assertEqual(len(_tags(base)), len(_tags(mirror)))
                for i, (b, m) in enumerate(zip(_tags(base), _tags(mirror))):
                    with self.subTest(tag=i):
                        self.assertNotEqual(b, m, "every stereocentre must invert")

    def test_achiral_fixture_mirror_is_IDENTICAL(self):
        """The other direction: an achiral molecule's mirror MUST encode identically.

        This is the over-sensitivity guard, and it is the one a naive "mirrors must always
        differ" rule gets backwards. ROGYAO_comp_0 is achiral by the torsion-aware oracle
        (RIGID_ACHIRAL, mirror superimposes at 0.423 A over 4 automorphisms), so an encoder
        that distinguished it would be manufacturing stereochemistry that is not there --
        a false positive, which the round trip would then fail on forever.
        """
        path = _fixture(_UNSTABLE)
        base = _encode(path, stable=True)
        with tempfile.TemporaryDirectory() as tmp:
            mirror = _encode(_mirrored(path, tmp), stable=True)
        self.assertEqual(
            base,
            mirror,
            f"{_UNSTABLE} is achiral (oracle: RIGID_ACHIRAL, 0.423 A) so its mirror must "
            "encode identically; a difference here means the encoder invented chirality",
        )

    def test_mirror_of_an_rr_fixture_is_ss(self):
        for name in _RR_FIXTURES:
            with self.subTest(fixture=name):
                with tempfile.TemporaryDirectory() as tmp:
                    mirror = _mirrored(_fixture(name), tmp)
                    self.assertEqual(_cip_labels(_encode(mirror, stable=True)), ["S", "S"])

    def test_tag_count_is_preserved(self):
        """The lever must never make the string stable by emitting fewer tags."""
        for name in (_UNSTABLE,) + _RR_FIXTURES:
            with self.subTest(fixture=name):
                path = _fixture(name)
                self.assertEqual(
                    len(_tags(_encode(path, stable=True))),
                    len(_tags(_encode(path, stable=False))),
                )


class TestLeverIsOffByDefault(unittest.TestCase):
    def test_unset_env_is_byte_identical_to_explicit_off(self):
        path = _fixture(_UNSTABLE)
        env = dict(os.environ)
        env.pop(_LEVER, None)
        with mock.patch.dict(os.environ, env, clear=True):
            unset = XYZToSMILES().convert(path)
        self.assertEqual(unset, _encode(path, stable=False))

    def test_lever_does_not_change_the_ORIGINAL_ordering(self):
        """On the input's own atom order the lever should be a no-op.

        This replaces an assertion that required ON != OFF here, which had the lever's
        contract backwards. The lever re-derives a chiral tag whose parity was recorded
        against a neighbour ordering the fragment rebuild destroyed -- so it changes the
        answer for *renumbered* presentations, and on the original ordering there is
        usually nothing to correct. ON == OFF here is the desirable property: no
        gratuitous change to the string everyone already has.

        Non-inertness is proved where it actually lives, by
        test_lever_off_reproduces_the_defect plus
        test_lever_on_is_byte_stable_under_renumbering.
        """
        path = _fixture(_UNSTABLE)
        self.assertEqual(
            _encode(path, stable=False),
            _encode(path, stable=True),
            "the lever changed the encoding of the ORIGINAL atom ordering; it is meant to "
            "repair renumbered presentations, not to rewrite the canonical answer",
        )


if __name__ == "__main__":
    unittest.main()
