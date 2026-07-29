"""``OIN_RESONANCE_DONOR_FOLD`` -- the frozen-resonance widening of the donor fold (v0.4.14).

WHAT IS BEING PINNED, AND WHY EACH TEST EXISTS
==============================================
The lever widens ``_donor_swap_permutations``' condition (a) so that two donors distinguished
only by which one the perceiver wrote as the ketone are exchangeable. Three properties have to
hold, and only the first is what the lever is *for*:

1. it merges genuine resonance pairs (acac, carboxylate, sulfonate);
2. it still REFUSES constitutional differences (ester, ether/ketone, different elements) --
   a widening that merged those would be folding two chemically distinct donors together;
3. it is byte-identical to v0.4.13 when off.

Property 3 is the one that would fail silently. ``_merge_classes`` is a union-find over two key
functions, and the second returns ``None`` with the lever off; if that ever stopped being true
the partition could change on the DEFAULT path while every one of this project's gates stayed
green, because the emitted difference is a slot index that the comparison key folds.

⚠ These are FIXTURE tests and their coverage is exactly three ligand motifs. They cannot tell
you the lever is safe on the corpus -- ``tools/mirror_audit_donor_fold.py`` is the instrument
for that, and v0.4.13 recorded ARM 1 passing byte-identically across a promotion because 0 of
its 62 fixtures were movers. State coverage, not verdicts.
"""

import os
import unittest
from unittest import mock

from oinsmiles.oin.canonical_slots import _skeleton_ranks, canonicalize_oin_slots
from oinsmiles.oin.compare import _parse_fragment
from oinsmiles.oin.levers import default_on, held_off, lever_enabled

LEVER = "OIN_RESONANCE_DONOR_FOLD"
FOLD = "OIN_CANONICAL_DONOR_FOLD"
VETO = "OIN_FOLD_PARITY_VETO"


def _ranks_merge(smiles, idxs, skeleton):
    """Are the atoms at ``idxs`` in one class of ``smiles``' fragment (or its skeleton)?"""
    from rdkit import Chem

    mol = _parse_fragment(smiles)
    assert mol is not None, f"fixture {smiles!r} does not parse -- fix the fixture, not the test"
    ranks = _skeleton_ranks(mol) if skeleton else list(Chem.CanonicalRankAtoms(mol, False))
    assert ranks is not None, f"skeleton ranking failed on {smiles!r}"
    return len({ranks[i] for i in idxs}) == 1


class TestSkeletonMergesResonanceOnly(unittest.TestCase):
    """The equivalence notion itself, independent of any OIN string."""

    #: ``(label, fragment SMILES, donor atom indices)``. Chosen so each pair is the SAME
    #: element -- the interesting question is never "does it tell O from N".
    MERGE = [
        ("acac (ketone/enol)", "CC(=O)C=C(C)O", (2, 6)),
        ("carboxylate", "OC(=O)C", (0, 2)),
        ("sulfonate", "OS(=O)(=O)c1ccccc1", (0, 2)),
    ]
    REFUSE = [
        ("ester -O- vs =O", "CCOC(C)=O", (2, 5)),
        ("ether vs ketone", "COCCC=O", (1, 5)),
        ("amide N vs O", "NC(C)=O", (0, 3)),
    ]

    def test_resonance_pairs_merge_under_the_skeleton(self):
        for label, smi, idxs in self.MERGE:
            with self.subTest(label):
                self.assertFalse(
                    _ranks_merge(smi, idxs, skeleton=False),
                    f"{label}: the STRICT ranking should still separate these -- if it does not, "
                    "this fixture no longer demonstrates the defect the lever exists for",
                )
                self.assertTrue(
                    _ranks_merge(smi, idxs, skeleton=True),
                    f"{label}: skeleton ranking must merge a resonance pair",
                )

    def test_constitutional_differences_are_still_refused(self):
        for label, smi, idxs in self.REFUSE:
            with self.subTest(label):
                self.assertFalse(
                    _ranks_merge(smi, idxs, skeleton=True),
                    f"{label}: the skeleton erases bond-order bookkeeping ONLY. Merging these "
                    "would fold two chemically distinct donors onto one slot class.",
                )

    def test_the_skeleton_keeps_chiral_tags(self):
        """Structural assertion: flattening erases bond orders, charges and H counts -- not tags.

        Asserted on the skeleton molecule directly rather than through a ranking, because a
        ranking result can agree with this claim for unrelated reasons.
        """
        from rdkit import Chem

        mol = _parse_fragment("OC[C@H](C)[C@@H](C)CO")
        before = [str(a.GetChiralTag()) for a in mol.GetAtoms()]
        rw = Chem.RWMol(mol)
        for b in rw.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in rw.GetAtoms():
            a.SetFormalCharge(0)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(0)
            a.SetIsAromatic(False)
        skel = rw.GetMol()
        Chem.SanitizeMol(skel, Chem.SANITIZE_SYMMRINGS | Chem.SANITIZE_ADJUSTHS)
        self.assertEqual(before, [str(a.GetChiralTag()) for a in skel.GetAtoms()])

    def test_the_skeleton_ranking_consumes_those_tags(self):
        """...and the ranking is actually SENSITIVE to them, not merely carrying them along.

        Two fragments differing in exactly one chiral tag must not produce the same rank vector.
        Without this, ``test_the_skeleton_keeps_chiral_tags`` would be satisfied by tags that
        ``CanonicalRankAtoms`` ignores.

        ⚠ WHAT THIS DOES **NOT** ESTABLISH, stated because the tempting reading is wrong.
        It does not mean the widening cannot merge two stereochemically distinguishable donors.
        Measured on these two fixtures: the C2-symmetric ``(R,R)`` diol does **not** merge its two
        arms (over-conservative -- a missed fold, the safe direction), while the meso ``(R,S)``
        one **does** -- its arms are enantiotopic, and folding them is exactly the reflection
        ``fold_parity``'s veto exists to police. The guard against that case is the veto, not this
        ranking. What this ranking supplies is only that the widening does not DISCARD the
        chirality information the v0.4.11 strict ranking already used.
        """
        rr = _skeleton_ranks(_parse_fragment("OC[C@H](C)[C@H](C)CO"))
        rs = _skeleton_ranks(_parse_fragment("OC[C@H](C)[C@@H](C)CO"))
        self.assertNotEqual(rr, rs, "the skeleton ranking must consume chiral tags")


class TestLeverMovesTheStringOnlyWhenOn(unittest.TestCase):
    """End-to-end on the two motifs the corpus measurement named, through the real post-pass."""

    #: Real ``key_equal/slot_renumber`` pairs from ``results-v0.4.8-honest`` -- the two strings
    #: are the same complex encoded from the input XYZ and from its round trip.
    PAIRS = {
        "ALAJON_comp_0 (acac)": (
            "[Pt_SPL].CC(=O{0})C=C(C)O{1}.Cc1cc(-c2c{2}ccc3ccccc23)n{3}c2ccccc12",
            "[Pt_SPL].CC(=O{1})C=C(C)O{0}.Cc1cc(-c2c{2}ccc3ccccc23)n{3}c2ccccc12",
        ),
        "AROKUP_comp_0 (sulfonate)": (
            "[Ru_OCT].O{0}S(=O)(=O{2})c1ccccc1P{4}(C1CCCCC1)C1CCCCC1."
            "S{1}c1c(-c2ccccc2)cc(-c2ccccc2)cc1-c1ccccc1.[CH]{3}c1ccccc1O{5}C(C)C",
            "[Ru_OCT].O{2}S(=O)(=O{0})c1ccccc1P{4}(C1CCCCC1)C1CCCCC1."
            "S{1}c1c(-c2ccccc2)cc(-c2ccccc2)cc1-c1ccccc1.[CH]{3}c1ccccc1O{5}C(C)C",
        ),
    }

    def _canon(self, s, resonance):
        env = {FOLD: "1", LEVER: "1" if resonance else "0"}
        with mock.patch.dict(os.environ, env):
            return canonicalize_oin_slots(s)

    def test_off_leaves_the_pair_divergent(self):
        for label, (a, b) in self.PAIRS.items():
            with self.subTest(label):
                self.assertNotEqual(
                    self._canon(a, False),
                    self._canon(b, False),
                    "with the lever OFF these must still differ, or the fixture is stale and "
                    "this test is measuring nothing",
                )

    def test_on_converges_the_pair(self):
        for label, (a, b) in self.PAIRS.items():
            with self.subTest(label):
                self.assertEqual(self._canon(a, True), self._canon(b, True))

    def test_idempotent(self):
        for label, (a, _b) in self.PAIRS.items():
            with self.subTest(label):
                once = self._canon(a, True)
                self.assertEqual(once, self._canon(once, True))

    def test_unset_is_byte_identical_to_the_lever_being_ON(self):
        """Property 3, restated for the v0.4.14 promotion.

        Before promotion this asserted ``off == unset``. It now asserts ``unset == on``, and the
        inversion is the point: leaving the lever alone must give the SHIPPED answer, at the
        string level and not merely in the registry. ``TestLeverRegistration`` checks
        ``default_on()``; that is a statement about a frozenset, and a lever can be listed there
        while a call site still spells its own default and reaches the old code path.
        """
        for label, (a, b) in self.PAIRS.items():
            for s in (a, b):
                with self.subTest(f"{label} {s[:24]}"):
                    with mock.patch.dict(os.environ, {FOLD: "1", LEVER: "1"}):
                        on = canonicalize_oin_slots(s)
                    env = {k: v for k, v in os.environ.items() if k != LEVER}
                    env[FOLD] = "1"
                    with mock.patch.dict(os.environ, env, clear=True):
                        unset = canonicalize_oin_slots(s)
                    self.assertEqual(on, unset)

    def test_explicit_off_still_reaches_the_v0413_answer(self):
        """The opt-out has to work, and ``"0"`` is the spelling that historically did not.

        ``os.environ.get(name)`` returns ``"0"`` as a truthy string, so a call site using it
        would ENABLE the lever here and this test would fail -- which is the only cheap way to
        catch that, since the wrong spelling is silent at runtime.
        """
        for label, (a, b) in self.PAIRS.items():
            with self.subTest(label):
                with mock.patch.dict(os.environ, {FOLD: "1", LEVER: "0"}):
                    self.assertNotEqual(canonicalize_oin_slots(a), canonicalize_oin_slots(b))

    def test_cannot_fire_with_the_donor_fold_off(self):
        """It widens a candidate set that does not exist unless the fold is on."""
        for label, (a, b) in self.PAIRS.items():
            with self.subTest(label):
                with mock.patch.dict(os.environ, {FOLD: "0", LEVER: "1"}):
                    self.assertNotEqual(canonicalize_oin_slots(a), canonicalize_oin_slots(b))


class TestLeverRegistration(unittest.TestCase):
    def test_the_shipped_default_is_on(self):
        """Pins the v0.4.14 promotion itself, so a silent revert is a test failure."""
        self.assertIn(LEVER, default_on())
        env = {k: v for k, v in os.environ.items() if k != LEVER}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertTrue(lever_enabled(LEVER), f"{LEVER} ships ON since v0.4.14")

    def test_is_no_longer_listed_as_held_off(self):
        self.assertNotIn(
            LEVER,
            held_off(),
            "_HELD_OFF is for levers that are deliberately opt-in. A promoted lever listed "
            "there makes held_off() lie about the shipped configuration.",
        )

    def test_read_through_lever_enabled_never_os_environ_get(self):
        """``os.environ.get(name)`` returns ``"0"``, which is truthy and ENABLES the lever.

        This project has paid 23 test failures across two promotions for that. Pinned by
        source inspection because the wrong spelling is silent at runtime.
        """
        import pathlib

        src = pathlib.Path(__file__).resolve().parents[2] / "src" / "oinsmiles"
        offenders = [
            str(p.relative_to(src))
            for p in src.rglob("*.py")
            if f'os.environ.get("{LEVER}"' in p.read_text()
        ]
        self.assertEqual(offenders, [], f"{LEVER} must be read via lever_enabled()")


class TestResonanceFoldInheritsTheVetoCoupling(unittest.TestCase):
    """The widening is subject to the SAME coupling invariant as the fold it widens.

    ``OIN_CANONICAL_DONOR_FOLD`` ON with ``OIN_FOLD_PARITY_VETO`` OFF collapses enantiomers, and
    no gate arm, golden or comparison key can see it. This lever makes that candidate set
    *larger*, so promoting it while the veto is off would enlarge exactly the damage v0.4.11
    refuted. It can never be promoted ahead of the veto.
    """

    def test_cannot_ship_on_while_the_veto_ships_off(self):
        on = default_on()
        if LEVER in on:
            self.assertIn(
                VETO,
                on,
                f"{LEVER} widens the donor fold's candidate set. Promoting it with "
                f"{VETO} off enlarges the enantiomer collapse v0.4.11 refuted.",
            )
            self.assertIn(FOLD, on, f"{LEVER} is a no-op without {FOLD}")


if __name__ == "__main__":
    unittest.main()
