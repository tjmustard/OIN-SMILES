"""The ``Atom count mismatch`` gate is load-bearing, and this pins why.

THE QUESTION THIS SETTLES
=========================
The harness fails a round trip when the generated XYZ has a different atom count from the
input. At N = 5000 that fires on **18 molecules**, and every one of them *already passed the
string comparison* -- the gate runs after it (``test_dataset_roundtrip.py``: key comparison,
then RMSD, then atom count). So the reasonable suspicion was that the gate is a third error
direction: 18 correct round trips scored as failures.

It is not. Measured over all 18 (``tools/atom_count_provenance.py``):

* the element delta is **hydrogen only** in 18/18, always ``+1`` or ``+2``, always gained;
* **8 of the 18 re-encode -- independently, from coordinates alone -- to the SAME OIN as the
  input.**

That second figure is the finding. For those 8, a structure carrying two extra hydrogens and
the original structure produce byte-identical OIN. No string comparison can separate them:
not the scored one, not the honest one, not the key. The notation is not injective over
hydrogen count for this class, and the atom-count gate is the only instrument in the harness
that can see it.

The contract this project exists to satisfy is **lossless** 3D <-> 1D. A round trip that
invents two hydrogens has broken it. Failing those molecules is correct.

⚠ WHY ``indep_key_match`` CANNOT BE USED AS EVIDENCE HERE. It is tempting to read "the key
matches, so the structures agree". The key is lossy -- it strips metal ``@``, folds slot drift
and folds the axial token -- and this lane establishes it folds hydrogen count too. A lossy key
must never be reused as an acceptance predicate for an axis it folds, which is exactly the
mistake this test exists to prevent someone repeating.

Already refuted, do not re-run: ``OIN_H_FAITHFUL`` was built for this class and bought nothing
(A/B over the 45-molecule population: match 8 / mismatch 37 with the lever off, and the
identical 8 / 37 with it on). See ``docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md``.
"""

import collections
import math
import os
import unittest

from oinsmiles import XYZToSMILES
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.oin.coordination import coordination_report, parse_xyz

_FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "atom_count_gate")

#: ``smiles_1`` as the v0.4.6 sweep recorded it. Hard-coded so the assertion fails if the
#: encoder moves -- an expectation recomputed from the encoder under test cannot detect that.
_XAKCAP_OIN = "[Mo_OCT].CN(C)CS{0}.S{1}=c1ccccn{3}1.C{2}#O.C{4}#O.c1ccc(P{5}(c2ccccc2)c2ccccc2)cc1"


def _read(name, kind):
    with open(os.path.join(_FIX, f"{name}_{kind}.xyz")) as fh:
        return fh.read()


def _counts(text):
    symbols, _ = parse_xyz(text)
    return collections.Counter(symbols)


class TestAtomCountGateIsLoadBearing(unittest.TestCase):
    """``XAKCAP_comp_0`` -- the sharpest case: invisible to every other instrument."""

    NAME = "XAKCAP_comp_0"

    def test_the_generated_structure_has_two_extra_hydrogens(self):
        delta = _counts(_read(self.NAME, "generated"))
        delta.subtract(_counts(_read(self.NAME, "input")))
        self.assertEqual(
            {k: v for k, v in delta.items() if v},
            {"H": 2},
            "the discrepancy must be hydrogen and nothing else",
        )

    def test_the_oin_string_cannot_see_it(self):
        """The decisive assertion. Two structures, 61 and 63 atoms, one OIN string."""
        honest = XYZToSMILES().convert(os.path.join(_FIX, f"{self.NAME}_generated.xyz"))
        self.assertEqual(
            _XAKCAP_OIN,
            honest,
            "a 63-atom structure must still encode to the 61-atom input's OIN -- if this "
            "ever stops being true the notation gained hydrogen fidelity and this lane's "
            "verdict needs re-deriving",
        )
        self.assertEqual(
            canonical_roundtrip_key(_XAKCAP_OIN),
            canonical_roundtrip_key(honest),
            "the key folds it too -- so `indep_key_match` is not evidence of agreement",
        )

    def test_the_coordination_probe_cannot_see_it_either(self):
        """Not a metal-coordination defect: the sphere is intact. Third instrument, blind."""
        report = coordination_report(_read(self.NAME, "input"), _read(self.NAME, "generated"))
        self.assertIsNot(report["intact"], False)

    def test_the_extra_hydrogens_are_not_a_dissociated_fragment(self):
        """Both added H sit on a heavy atom, so this is protonation and not stray atoms.

        Worth pinning separately: "two extra H" would be a much less interesting finding if
        they were floating unbonded, and a reader is entitled to know which it is.
        """
        symbols, coords = parse_xyz(_read(self.NAME, "generated"))
        heavy = [i for i, s in enumerate(symbols) if s != "H"]
        for i, s in enumerate(symbols):
            if s != "H":
                continue
            nearest = min(math.dist(coords[i], coords[j]) for j in heavy)
            self.assertLess(nearest, 1.3, f"hydrogen {i} is not bonded to anything")


class TestAtomCountGatePopulationShape(unittest.TestCase):
    """``KIKSAB_comp_0`` -- the other half: the honest arm DOES catch this one.

    The 18 split 8 / 10 on whether independent re-perception notices. Pinning one from each
    side keeps the lane's conclusion from being read as "the gate catches everything" or as
    "the honest metric makes the gate redundant". Neither is true.
    """

    NAME = "KIKSAB_comp_0"

    def test_hydrogen_only_discrepancy(self):
        delta = _counts(_read(self.NAME, "generated"))
        delta.subtract(_counts(_read(self.NAME, "input")))
        self.assertEqual({k: v for k, v in delta.items() if v}, {"H": 2})

    def test_the_honest_arm_catches_this_one(self):
        scored = (
            "[Ru_OCT].CC(=O{0})C=C(C)O{2}.[CH]1[CH][CH]C(N{4}N=c2cccc(=NN{5}C3[CH][CH][CH]"
            "[CH][CH]3)n{1}2)[CH][CH]1.[Cl]{3}"
        )
        honest = XYZToSMILES().convert(os.path.join(_FIX, f"{self.NAME}_generated.xyz"))
        self.assertNotEqual(scored, honest, "independent re-perception is expected to diverge here")


if __name__ == "__main__":
    unittest.main()
