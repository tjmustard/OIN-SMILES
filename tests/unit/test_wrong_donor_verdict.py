"""The `DISTINCT_donors` verdict was a false positive — pin why (v0.4.5 Lane 9).

Lane 2's `tools/slot_drift_mechanism.py` classified 7 residual `slot_renumber` pairs as
`DISTINCT_donors`, i.e. "one of the two strings is WRONG about which atom is bound where".
Lane 9 settled all 7 from the 3D coordinates (`tools/wrong_donor_groundtruth.py`) and found
**0 soundness defects**: 4 had a bit-identical donor->vertex map and the other 3 differed by a
proper rotation of the coordination polyhedron, `|delta rssd| <= 1.2e-14`. See
`docs/WRONG_DONOR_v0.4.5.md`.

These tests use the **real emitted strings** from the actual 7, so they need no dataset (the
dataset is gitignored and absent from a worktree). They pin the two classifier mechanisms that
were repaired and the group-theoretic fact that makes a per-fragment test misleading in
principle, so neither the bugs nor the unsupportable soundness claim can come back.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "tools")))

from slot_drift_mechanism import (  # noqa: E402  (tools/ is not a package)
    _pair_fragments,
    _slot_to_atoms,
    atom_verdict,
    mechanism,
)

from oinsmiles.oin.canonical_slots import geometry_rotation_group

# --- ZACFER_comp_0: two equivalent Cp* rings bridged by Si-Si, on a TET Ti center. -----
# The rings swap slot integers. They are interchangeable by the ligand's own C2, so this is
# an automorphism -- but the old code keyed each slot on the marker's FIRST occurrence, which
# compared a methyl-bearing ring carbon against a silyl-bearing one and called them distinct.
ZACFER_BASE = (
    "[Ti_TET].Cc{1>}1c{1}(C)c{1}(C)c{1}([Si](C)(C)[Si](C)(C)c{0}2c{0}(C)c{0>}(C)c{0}(C)c{0}2C)"
    "c{1}1C.O{2}S(=O)(=O)C(F)(F)F.[OH]{3}"
)
ZACFER_GOT = (
    "[Ti_TET].Cc{0>}1c{0}(C)c{0}(C)c{0}([Si](C)(C)[Si](C)(C)c{1}2c{1}(C)c{1>}(C)c{1}(C)c{1}2C)"
    "c{0}1C.O{2}S(=O)(=O)C(F)(F)F.[OH]{3}"
)

# --- RUBTIS_comp_0: COD + a pyridyl-imine chelate on a SPL Rh center. -----------------
# The chelate's donors ARE inequivalent (imine N vs pyridyl n), so per-fragment the swap looks
# illegitimate. Measured from the 3D: the physical donor -> vertex map is IDENTICAL between the
# two encodings; the COD's two equivalent arms swapped as well, and the composite relabeling is
# (0 1)(2 3), a proper rotation of the square.
RUBTIS_BASE = (
    "[Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1.CC(=N{2}c1c(C(C)C)cccc1C(C)C)c1ccccn{3}1"
)
RUBTIS_GOT = (
    "[Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1.CC(=N{3}c1c(C(C)C)cccc1C(C)C)c1ccccn{2}1"
)


class TestBothStringsAreTheSameMolecule(unittest.TestCase):
    """These pairs are re-presentations of ONE structure, so any real defect is downstream."""

    def test_both_pairs_are_same_vcolor_identical(self):
        """First-stage class must stay `same_vcolor_identical` -- the premise of stage two."""
        self.assertEqual(mechanism(ZACFER_BASE, ZACFER_GOT), "same_vcolor_identical")
        self.assertEqual(mechanism(RUBTIS_BASE, RUBTIS_GOT), "same_vcolor_identical")


class TestHapticFirstOccurrenceFalsePositive(unittest.TestCase):
    """ZACFER: an eta ring's slot must be identified by ALL its atoms, not the first one."""

    def test_zacfer_is_an_automorphism(self):
        self.assertEqual(atom_verdict(ZACFER_BASE, ZACFER_GOT), "automorphism")

    def test_slot_markers_collect_every_occurrence(self):
        """A 5-atom eta group must contribute 5 atoms to its slot, not 1."""
        frag = ZACFER_BASE.split(".")[1]
        by_slot = _slot_to_atoms(frag)
        self.assertEqual(len(by_slot[0]), 5, "eta5 ring should carry 5 atoms at its slot")
        self.assertEqual(len(by_slot[1]), 5, "eta5 ring should carry 5 atoms at its slot")


class TestFragmentsPairedByBodyNotPosition(unittest.TestCase):
    """The post-pass re-derives fragment order FROM the slots, so position is not stable."""

    def test_reordered_fragments_still_pair_on_body(self):
        base = ["[Cl]{0}", "N{1}CC", "N{2}CC"]
        got = ["N{2}CC", "[Cl]{0}", "N{1}CC"]  # same bodies, different order
        pairs = _pair_fragments(base, got)
        self.assertIsNotNone(pairs)
        for a, b in pairs:
            self.assertEqual(
                a.replace("{0}", "").replace("{1}", "").replace("{2}", ""),
                b.replace("{0}", "").replace("{1}", "").replace("{2}", ""),
                "fragments must be paired by body text",
            )

    def test_absent_body_is_reported_not_mispaired(self):
        self.assertIsNone(_pair_fragments(["N{0}CC"], ["O{0}CC"]))


class TestVerdictMakesNoSoundnessClaim(unittest.TestCase):
    """`DISTINCT_donors` asserted a wrong string. Measured 0/7. The name must not come back."""

    def test_rubtis_verdict_is_the_local_heuristic_name(self):
        verdict = atom_verdict(RUBTIS_BASE, RUBTIS_GOT)
        self.assertEqual(verdict, "distinct_donors_LOCAL")
        self.assertNotIn("DISTINCT_donors", verdict)

    def test_module_does_not_reintroduce_the_soundness_claim(self):
        import slot_drift_mechanism as sdm

        with open(sdm.__file__) as fh:
            src = fh.read()
        self.assertNotIn(
            'return "DISTINCT_donors"',
            src,
            "the verdict cannot claim a string is wrong -- see docs/WRONG_DONOR_v0.4.5.md",
        )


class TestWhyAPerFragmentTestIsMisleading(unittest.TestCase):
    """The group-theoretic fact behind RUBTIS, pinned so the reasoning survives.

    A per-fragment test sees only `(2 3)` -- swapping the chelate's two inequivalent donors --
    which is NOT a rotation of the square, so it flags the pair. The relabeling that actually
    happened is `(0 1)(2 3)`: the COD's two equivalent arms swapped too, and *that* IS a proper
    rotation. Composition across fragments is invisible to a fragment-local question.
    """

    def test_chelate_swap_alone_is_not_a_rotation_but_the_composite_is(self):
        group = {tuple(g) for g in geometry_rotation_group("SPL")}
        chelate_only = (0, 1, 3, 2)  # (2 3)
        composite = (1, 0, 3, 2)  # (0 1)(2 3)
        self.assertNotIn(chelate_only, group, "(2 3) alone must not be a proper rotation of SPL")
        self.assertIn(composite, group, "(0 1)(2 3) must be a proper rotation of SPL")

    def test_every_rotation_is_proper_so_chirality_cannot_be_folded(self):
        """Folding over this group can never merge enantiomers -- the Y2 axial hazard."""
        import numpy as np

        from oinsmiles.oin.canonical_slots import GEOMETRY_VERTICES

        for geo in ("OCT", "TET", "SPL", "PBP"):
            verts = np.asarray(GEOMETRY_VERTICES[geo], dtype=float)
            verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)
            for perm in geometry_rotation_group(geo):
                # The Gram matrix is preserved by construction; assert the realizing map is a
                # proper rotation for the spanning (rank-3) templates.
                if np.linalg.matrix_rank(verts, tol=1e-6) < 3:
                    continue
                basis: list[int] = []
                for i in range(len(verts)):
                    if np.linalg.matrix_rank(verts[basis + [i]], tol=1e-6) == len(basis) + 1:
                        basis.append(i)
                    if len(basis) == 3:
                        break
                rot = verts[[perm[b] for b in basis]].T @ np.linalg.inv(verts[basis].T)
                self.assertGreater(
                    np.linalg.det(rot), 0, f"{geo} group must contain proper rotations only"
                )


if __name__ == "__main__":
    unittest.main()
