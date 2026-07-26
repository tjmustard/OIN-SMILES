"""Coordination-integrity diagnostic — `oin/coordination.py`.

This exists because the round-trip metric credits structures whose ligands have come off: it
scores through `gen_result.mol`, the generator's own bond graph. Measured false-positive rate on
the default path is 9.6% overall and 28.1% on haptic molecules
(`docs/METRIC_FALSE_POSITIVES.md`).

Synthetic geometries rather than fixtures, so each test pins exactly one behaviour and the
distances are visible in the test itself.
"""

import unittest

from oinsmiles.oin.coordination import (
    MARGINAL_BAND,
    coordination_report,
    metal_contacts,
    metal_indices,
    parse_xyz,
)


def xyz(atoms) -> str:
    """`[(symbol, x, y, z), ...]` -> XYZ file content."""
    lines = [str(len(atoms)), "test"]
    lines += [f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, x, y, z in atoms]
    return "\n".join(lines) + "\n"


def ferrocene_like(m_c: float) -> str:
    """Fe with ten carbons at radius `m_c`, five above and five below — a sandwich.

    Real ferrocene Fe-C is ~2.05 A. The Fe/C contact cutoff is ~2.53 A, so `m_c=2.05` is bonded
    and `m_c=2.90` (what the generator produced for FIYHUT_comp_0) is not.
    """
    import math

    atoms = [("Fe", 0.0, 0.0, 0.0)]
    for sign in (1, -1):
        for k in range(5):
            a = 2 * math.pi * k / 5
            # place on a cone so the straight-line Fe-C distance is exactly m_c
            r = m_c * 0.8
            z = sign * (m_c**2 - r**2) ** 0.5
            atoms.append(("C", r * math.cos(a), r * math.sin(a), z))
    return xyz(atoms)


class TestParsing(unittest.TestCase):
    def test_roundtrips_a_simple_file(self):
        syms, coords = parse_xyz(xyz([("Fe", 0, 0, 0), ("C", 2.0, 0, 0)]))
        self.assertEqual(syms, ["Fe", "C"])
        self.assertAlmostEqual(coords[1][0], 2.0)

    def test_malformed_yields_empty_rather_than_raising(self):
        for bad in ("", "not-a-number\ncomment\n", "2\ncomment\nFe 0 0 0\n"):
            with self.subTest(bad=bad[:14]):
                syms, _ = parse_xyz(bad)
                self.assertEqual(syms, [])

    def test_finds_the_metal_and_ignores_ligand_atoms(self):
        syms, _ = parse_xyz(xyz([("C", 0, 0, 0), ("Fe", 2, 0, 0), ("N", 4, 0, 0)]))
        self.assertEqual(metal_indices(syms), [1])


class TestContacts(unittest.TestCase):
    def test_counts_only_what_is_inside_the_cutoff(self):
        syms, coords = parse_xyz(xyz([("Fe", 0, 0, 0), ("C", 2.0, 0, 0), ("C", 5.0, 0, 0)]))
        counts, contacts = metal_contacts(syms, coords, 0)
        self.assertEqual(counts, {"C": 1})
        self.assertEqual(len(contacts), 1)

    def test_contacts_are_sorted_by_distance(self):
        syms, coords = parse_xyz(xyz([("Fe", 0, 0, 0), ("C", 2.3, 0, 0), ("N", 1.9, 0, 0)]))
        _counts, contacts = metal_contacts(syms, coords, 0)
        self.assertEqual([c[1] for c in contacts], ["N", "C"])


class TestVerdict(unittest.TestCase):
    def test_identical_geometry_is_intact(self):
        s = ferrocene_like(2.05)
        rep = coordination_report(s, s)
        self.assertTrue(rep["intact"])
        self.assertEqual(rep["metals"][0]["n_contacts_gen"], 10)

    def test_detached_sandwich_is_flagged(self):
        """The FIYHUT_comp_0 failure: both rings ~0.85 A too far, 10 bonded carbons -> 0."""
        rep = coordination_report(ferrocene_like(2.05), ferrocene_like(2.90))
        self.assertIs(rep["intact"], False)
        m = rep["metals"][0]
        self.assertEqual(m["n_contacts_in"], 10)
        self.assertEqual(m["n_contacts_gen"], 0)
        self.assertEqual(m["lost"], {"C": 10})
        self.assertEqual(m["lost_beyond_band"], 10)
        self.assertIn("lost", rep["reason"])

    def test_a_loss_entirely_at_the_boundary_is_NOT_charged(self):
        """Contacts a hair outside the cutoff are inconclusive, not a degradation.

        Added on evidence: a raw loss verdict flagged 45 genuine passes on the 633-molecule
        validation set, 36 of them with a contact within MARGINAL_BAND. Distinguishing them cut
        false alarms 7.9% -> 3.7% for one point of recall.
        """
        # Fe/C cutoff ~2.53; nudge just past it by less than the band.
        syms, _ = parse_xyz(ferrocene_like(2.05))
        cutoff = metal_contacts(*parse_xyz(ferrocene_like(2.05)), 0)[1][0][3]
        rep = coordination_report(
            ferrocene_like(2.05), ferrocene_like(cutoff + MARGINAL_BAND * 0.4)
        )
        self.assertIsNot(rep["intact"], False, "a boundary drift must not read as degraded")
        self.assertTrue(rep["boundary_only"])
        del syms

    def test_unassessable_is_None_and_never_True(self):
        """A silent True is exactly how the metric this replaces fails."""
        for a, b in (
            ("", ferrocene_like(2.05)),
            (xyz([("C", 0, 0, 0), ("H", 1, 0, 0)]), xyz([("C", 0, 0, 0), ("H", 1, 0, 0)])),
        ):
            with self.subTest():
                rep = coordination_report(a, b)
                self.assertIsNone(rep["intact"])
                self.assertIsNotNone(rep["reason"])

    def test_gaining_a_contact_is_reported_but_not_charged_as_loss(self):
        """A gain is not reliably bad -- a genuine pass in the corpus gained 2 (Mo 6->8).

        Recorded so a caller can see over-coordination, without this module asserting a threshold
        it has no evidence for. 4 of 61 known false positives are gain-driven and are a documented
        scope limit, not a silent miss.
        """
        base = xyz([("Fe", 0, 0, 0), ("C", 2.0, 0, 0)])
        more = xyz([("Fe", 0, 0, 0), ("C", 2.0, 0, 0), ("C", 0, 2.0, 0)])
        rep = coordination_report(base, more)
        self.assertTrue(rep["intact"])
        self.assertEqual(rep["metals"][0]["gained"], {"C": 1})


if __name__ == "__main__":
    unittest.main()
