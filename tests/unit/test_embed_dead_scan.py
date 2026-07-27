"""Guards for v0.4.10's removal of the discarded ``.index()`` scan in ``get_embedding``.

WHAT WAS REMOVED, AND WHY IT WAS NOT A NO-OP
=============================================
``generator3d/embed.py``'s ``get_embedding`` opened its outer loop with

    for alternative_ace_mol in alternative_ace_mol_list:
        alternative_ace_mol_list.index(alternative_ace_mol)   # result DISCARDED

``list.index`` compares with ``Molecule.__eq__`` -> ``is_same_molecule`` ->
``get_c_eig_list`` -> ``numpy.linalg.eig``, so this eigendecomposed a Coulomb matrix per
candidate per outer iteration and threw the answer away. v0.4.9 costed one profile line of
it -- ``numpy.linalg.eig``, 15.99 s of a 72.02 s generation, 3711 index calls and 198
eigendecompositions (``docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md`` §5). The A/B on
the whole call measured **-50.2% on ``CAHQEJ_comp_0``**, because the scan also rebuilds a
Coulomb matrix on both operands of every comparison that survives the atom-count check.
It measured **nothing on ``FOSNEI_comp_0``** -- the cost is bimodal by molecule.

WHY THESE TESTS EXIST RATHER THAN JUST THE GATE
================================================
The byte-identity gate proves the removal changed nothing *today*. It cannot explain
*why*, so it cannot stop someone from re-introducing the dependency later. The deletion is
sound only while the whole comparison call graph is side-effect-free -- in particular
while ``get_c_eig_list`` keeps recomputing instead of populating ``self.c_eig_list``. It
already reads ``if self.c_eig_list is None:`` and then never assigns, so adding the
"obvious" one-line cache there is a live temptation. The moment someone does, the deleted
call becomes an observable cache-warmer and its removal stops being byte-identical.

These tests pin that property directly, so that change fails here with an explanation
rather than surfacing as an unattributable string diff in a 10-CPU-hour corpus run.
"""

import pathlib
import unittest

import numpy as np

from oinsmiles.generator3d.chem import Molecule

EMBED_PY = (
    pathlib.Path(__file__).resolve().parents[2] / "src" / "oinsmiles" / "generator3d" / "embed.py"
)


def _snapshot(mol):
    """A comparable snapshot of every attribute on *mol*.

    numpy arrays are compared by content, so a recomputed-but-equal matrix does not read
    as a mutation while an actual cache fill does.
    """
    out = {}
    for k, v in vars(mol).items():
        if isinstance(v, np.ndarray):
            out[k] = ("ndarray", v.tobytes(), v.shape)
        elif isinstance(v, list):
            out[k] = ("list", len(v))
        else:
            out[k] = ("plain", repr(v))
    return out


class TestComparisonCallGraphIsPure(unittest.TestCase):
    """The property that makes the deletion byte-identical."""

    def test_get_c_eig_list_does_not_populate_its_own_cache(self):
        """If this starts failing, the deleted ``.index()`` call was warming a cache."""
        mol = Molecule("CCO")
        self.assertIsNone(mol.c_eig_list)
        first = mol.get_c_eig_list()
        self.assertIsNone(
            mol.c_eig_list,
            "get_c_eig_list now populates self.c_eig_list. That makes the discarded "
            "list.index() scan removed in v0.4.10 an observable cache-warmer, so its "
            "removal is no longer byte-identical. Either revert the caching or re-derive "
            "the deletion's justification -- do not just update this assertion.",
        )
        second = mol.get_c_eig_list()
        np.testing.assert_array_equal(first, second)

    def test_is_same_molecule_mutates_neither_operand(self):
        a, b = Molecule("CCO"), Molecule("CCO")
        before_a, before_b = _snapshot(a), _snapshot(b)
        self.assertTrue(a.is_same_molecule(b, True))
        self.assertEqual(before_a, _snapshot(a))
        self.assertEqual(before_b, _snapshot(b))

    def test_eq_mutates_neither_operand(self):
        """``__eq__`` is what ``list.index`` actually calls."""
        a, b = Molecule("CCO"), Molecule("CCN")
        before_a, before_b = _snapshot(a), _snapshot(b)
        _ = a == b
        self.assertEqual(before_a, _snapshot(a))
        self.assertEqual(before_b, _snapshot(b))

    def test_index_over_a_member_cannot_raise(self):
        """The removed call could not have been a disguised guard.

        ``list.index`` short-circuits on identity, so searching for an element drawn from
        the list itself always finds it -- ValueError was unreachable and no control flow
        depended on it.
        """
        mols = [Molecule("CCO"), Molecule("CCN"), Molecule("CCC")]
        for m in mols:
            self.assertIsInstance(mols.index(m), int)


class TestDeadScanStaysRemoved(unittest.TestCase):
    """A source lint: cheap, and the only thing that stops a silent 50% regression."""

    def test_no_discarded_index_scan_in_embed(self):
        src = EMBED_PY.read_text()
        offenders = [
            (n, line.strip())
            for n, line in enumerate(src.splitlines(), 1)
            if line.strip() == "alternative_ace_mol_list.index(alternative_ace_mol)"
        ]
        self.assertEqual(
            offenders,
            [],
            "The discarded list.index() scan is back in get_embedding. Removing it "
            "measured -50.2% on CAHQEJ_comp_0 (3711 calls / 198 eigendecompositions), and "
            "its result is never used.",
        )


if __name__ == "__main__":
    unittest.main()
