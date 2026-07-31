"""Tests for the attachment guard on the RETURN path (v0.4.15 Lane 1, ``OIN_ATTACH_RETURN``).

v0.4.7 built a coordinate-only attachment predicate and wired it into ACCEPTANCE only.
``OIN_ATTACH_CHECK``'s own lever entry names the gap that left as its third residual class:
"GAVSED (acceptance rejected everything, and ``_select_by_geometry``'s fallback ranking never
consults this check -- the check guards ACCEPTANCE, not RETURN)". These tests pin the return-path
version, and above all they pin the two properties that make it safe to ship:

1. **Lever off is byte-identical.** ``_attach_rank`` returns 0 for every conformer, so every sort
   key and every scan collapses to its pre-lever order.
2. **Lever on never returns worse.** When no conformer holds its sites, the lowest-energy one --
   the pre-lever answer -- is still what comes back. Only ``OIN_ATTACH_RETURN_STRICT`` converts
   that into a failure, and that is a separate lever precisely because it CAN lower the headline.
"""

import os
import unittest

from rdkit import Chem

from oinsmiles.generation.metallogen_adapter import _attach_rank
from oinsmiles.oin.levers import held_off, lever_enabled


def _ptcl2n2(departed=False):
    """Square-planar Pt with 4 donors, all bonded. ``departed`` moves atom 4 out of range.

    The bonds are kept in BOTH cases on purpose: a detached ligand keeps its bond object, so a
    fixture that also removed the bond would not exercise the trap this predicate exists for.
    """
    rw = Chem.RWMol()
    for zn in (78, 17, 17, 7, 7):
        rw.AddAtom(Chem.Atom(zn))
    for j in (1, 2, 3, 4):
        rw.AddBond(0, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    last = (0.0, 6.0, 0.0) if departed else (0.0, -2.0, 0.0)
    for i, p in enumerate([(0, 0, 0), (2, 0, 0), (-2, 0, 0), (0, 2, 0), last]):
        conf.SetAtomPosition(i, tuple(float(v) for v in p))
    mol.AddConformer(conf)
    return mol


class _Lever:
    """Set/restore an env lever around a block, so a failure cannot leak into later tests."""

    def __init__(self, **levers):
        self.levers = levers
        self.prior = {}

    def __enter__(self):
        for k, v in self.levers.items():
            self.prior[k] = os.environ.get(k)
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return self

    def __exit__(self, *exc):
        for k, v in self.prior.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return False


class TestLeverRegistration(unittest.TestCase):
    def test_both_levers_are_registered_as_held_off_with_a_reason(self):
        reasons = held_off()
        for name in ("OIN_ATTACH_RETURN", "OIN_ATTACH_RETURN_STRICT"):
            self.assertIn(name, reasons, f"{name} must carry its rationale in the registry")

    def test_the_rationale_states_the_measured_target_and_the_measured_exposure(self):
        why = held_off()["OIN_ATTACH_RETURN"]
        self.assertIn("289", why, "the measured target belongs in the justification")
        self.assertIn("52", why, "the measured EXPOSURE belongs there too, not just the win")
        self.assertIn("POVPIA", why, "the known residual must be stated, not hidden")

    def test_both_levers_default_off(self):
        with _Lever(OIN_ATTACH_RETURN=None, OIN_ATTACH_RETURN_STRICT=None):
            self.assertFalse(lever_enabled("OIN_ATTACH_RETURN"))
            self.assertFalse(lever_enabled("OIN_ATTACH_RETURN_STRICT"))

    def test_a_zero_string_disables_rather_than_enables(self):
        """``os.environ.get`` returns ``"0"`` as truthy; that bug cost 23 test failures twice."""
        with _Lever(OIN_ATTACH_RETURN="0"):
            self.assertFalse(lever_enabled("OIN_ATTACH_RETURN"))


class TestAttachRankIsInertWhenOff(unittest.TestCase):
    """The byte-identity property, stated as a test rather than as a comment."""

    def test_a_detached_conformer_ranks_zero_with_the_lever_off(self):
        mol = _ptcl2n2(departed=True)
        with _Lever(OIN_ATTACH_RETURN=None):
            self.assertEqual(
                _attach_rank(mol),
                0,
                "off, every conformer must rank equal so the sort collapses to pre-lever order",
            )

    def test_an_attached_conformer_also_ranks_zero_with_the_lever_off(self):
        with _Lever(OIN_ATTACH_RETURN=None):
            self.assertEqual(_attach_rank(_ptcl2n2(departed=False)), 0)


class TestAttachRankDiscriminatesWhenOn(unittest.TestCase):
    def test_detached_ranks_worse_than_attached(self):
        with _Lever(OIN_ATTACH_RETURN="1"):
            attached = _attach_rank(_ptcl2n2(departed=False))
            detached = _attach_rank(_ptcl2n2(departed=True))
        self.assertEqual((attached, detached), (0, 1))
        self.assertLess(attached, detached, "the sort must prefer the attached conformer")

    def test_the_departed_ligand_still_has_its_bond(self):
        """The premise of the whole predicate: a ``GetBonds()`` test would pass this fixture."""
        mol = _ptcl2n2(departed=True)
        self.assertEqual(len(mol.GetAtomWithIdx(0).GetBonds()), 4)

    def test_nothing_to_judge_abstains_rather_than_demoting(self):
        """A conformer is demoted on evidence, never on ignorance -- errors stay PERMISSIVE."""
        with _Lever(OIN_ATTACH_RETURN="1"):
            self.assertEqual(_attach_rank(None), 0)

    def test_an_unevaluable_mol_abstains(self):
        """No conformer -> ``conformer_ligands_attached`` abstains True -> rank 0."""
        mol = Chem.RWMol()
        mol.AddAtom(Chem.Atom(78))
        with _Lever(OIN_ATTACH_RETURN="1"):
            self.assertEqual(_attach_rank(mol.GetMol()), 0)


class TestSortKeyOrdering(unittest.TestCase):
    """``scored`` is ``(attach_bad, clash, fit, rank)``; attachment must LEAD."""

    def test_an_attached_worse_fitting_conformer_outranks_a_detached_better_fitting_one(self):
        # (attach_bad, clash, fit, energy_rank) -- the detached one wins on every other term.
        detached_but_perfect = (1, 0, 0.0001, 0)
        attached_but_poor = (0, 0, 0.9, 7)
        self.assertLess(
            attached_but_poor,
            detached_but_perfect,
            "template fit means nothing if the ligand is not on the metal",
        )

    def test_with_the_lever_off_the_leading_term_cannot_reorder_anything(self):
        """Off, attach_bad is 0 everywhere, so the tuple order is decided by (clash, fit, rank)."""
        rows = [(0, 0, 0.9, 7), (0, 0, 0.0001, 0)]
        self.assertEqual(sorted(rows)[0], (0, 0, 0.0001, 0))


if __name__ == "__main__":
    unittest.main()
