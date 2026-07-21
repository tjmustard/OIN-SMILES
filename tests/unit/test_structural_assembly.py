"""Guards for A2 (v0.4.3 no-conformers wave): a complex that cannot be assembled
must surface a typed, diagnosable ``StructuralAssemblyError`` instead of being
silently absorbed into pool exhaustion (the generic ``no_conformers``).

The confirmed defect was the blanket ``except Exception`` in the embed attempt
loop of ``generate_3d_structures``. On the 7 real bucket-B molecules that fired
A0's ``pool.blanket_exception`` probe, the swallowed error was RDKit's
``AtomValenceException`` (a ``MolSanitizeException``): a dative donor drawn
covalently over-valences its atom -- e.g. FEJFAD_comp_0
``[Zn_TET].CN{0}(C)...N{1}(C)C...``, whose tertiary-amine N becomes 4-valent when
bonded to Zn ("Explicit valence for atom N, 4, is greater than permitted"). It
fired on all 250 attempts (deterministic), so the pool never filled and the
harness reported only ``MetalloGen failed to generate any conformers``. The
originally-hypothesised under-coordination path (``TypeError``/``IndexError`` in
embed's cmap construction) is covered too.

The loop body is monkeypatched, so these tests are deterministic and hermetic
(no real embedding, no g-xTB). Each ``surfaces_typed`` test fails against the
pre-fix code, which swallowed the error and returned ``[]``.
"""

import os
import sys
import unittest
from unittest import mock

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(os.path.join(_ROOT, "src"))
sys.path.append(os.path.join(_ROOT, "tools"))

import classify_failures  # noqa: E402
from rdkit import Chem  # noqa: E402

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles  # noqa: E402
from oinsmiles.generation.oin_parser import OINParser  # noqa: E402
from oinsmiles.generator3d import (  # noqa: E402
    StructuralAssemblyError,
    embed,
    generate_3d_structures,
)

# A trivial, always-parseable complex; the embed call is monkeypatched, so the
# only requirement is that om.get_om_from_modified_smiles() succeeds on it.
CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


def _raise_atom_valence(*a, **k):
    """Raise a genuine RDKit AtomValenceException (over-valent tertiary amine),
    the confirmed real cause on the 7 bucket-B molecules."""
    mol = Chem.MolFromSmiles("CN(C)(C)C", sanitize=False)  # 4-valent neutral N
    Chem.SanitizeMol(mol)  # raises AtomValenceException


class TestStructuralAssemblySurfaces(unittest.TestCase):
    def setUp(self):
        self.msmiles = _msmiles(CISPLATIN_OIN)

    def _run(self, side_effect):
        with mock.patch.object(embed, "get_embedding", side_effect=side_effect):
            return generate_3d_structures(self.msmiles, ff_params={"max_attempts": 5})

    def test_atom_valence_exception_surfaces_typed(self):
        """The confirmed real cause (an over-valent dative donor) surfaces typed."""
        with self.assertRaises(StructuralAssemblyError) as ctx:
            self._run(_raise_atom_valence)
        self.assertIn("structural error", str(ctx.exception))
        # The underlying RDKit cause is chained for the full traceback.
        self.assertIsInstance(ctx.exception.__cause__, Chem.AtomValenceException)

    def test_typeerror_surfaces_typed(self):
        """A None binding_site (TypeError) also surfaces as typed."""

        def structural_type_fail(*a, **k):
            raise TypeError("unsupported operand type(s) for -: 'NoneType' and 'int'")

        with self.assertRaises(StructuralAssemblyError):
            self._run(structural_type_fail)

    def test_indexerror_surfaces_typed(self):
        """An empty direction_vector (IndexError) also surfaces as typed."""

        def structural_index_fail(*a, **k):
            raise IndexError("list index out of range")

        with self.assertRaises(StructuralAssemblyError):
            self._run(structural_index_fail)

    def test_transient_exception_stays_no_conformers(self):
        """A genuinely transient (non-structural) failure keeps the old behavior:
        the pool exhausts to [] (generic no_conformers), NOT a typed error."""

        def transient_fail(*a, **k):
            raise RuntimeError("transient embed hiccup")

        self.assertEqual(
            self._run(transient_fail),
            [],
            "a transient (non-structural) failure must stay no_conformers, not become typed",
        )

    def test_none_return_stays_no_conformers(self):
        """An embed that just fails to validate (returns None) stays no_conformers."""
        self.assertEqual(self._run(lambda *a, **k: None), [])

    def test_mixed_structural_and_transient_stays_no_conformers(self):
        """If SOME attempts fail structurally but others only fail transiently
        (embed returned None), the complex is not *uniformly* unassemblable -- it
        must stay no_conformers ([]), never a typed error. Guards the real
        cisplatin regression: an over-valent bond-order *alternative* (option 1/2)
        must not mask a transient embed failure on a valid alternative (option 0)."""
        calls = {"n": 0}

        def mixed(*a, **k):
            calls["n"] += 1
            if calls["n"] % 2 == 1:
                _raise_atom_valence()  # odd calls: structural
            return None  # even calls: transient (embed ran, produced nothing)

        self.assertEqual(
            self._run(mixed),
            [],
            "a pool that is only partly structural must degrade to no_conformers",
        )

    def test_error_is_a_value_error(self):
        """Subclassing ValueError keeps existing `except ValueError` callers working."""
        self.assertTrue(issubclass(StructuralAssemblyError, ValueError))


class TestStructuralAssemblyClassified(unittest.TestCase):
    """The surfaced error must classify as a distinct reason, not no_conformers."""

    @staticmethod
    def _report(error):
        return {"status": "failed", "error": error, "smiles_1": CISPLATIN_OIN, "smiles_2": ""}

    def test_valence_cause_classifies_as_structural_assembly(self):
        err = (
            "Generation/Verification failed at UFF_1: StructuralAssemblyError: "
            "Could not assemble a valid 3D complex: every embed attempt failed with "
            "the same structural error (AtomValenceException: Explicit valence for "
            "atom # 2 N, 4, is greater than permitted)\n"
            "Traceback (most recent call last): ..."
        )
        cls, _ = classify_failures.classify(self._report(err))
        self.assertEqual(cls, "structural_assembly")

    def test_indexerror_cause_still_classifies_as_structural_assembly(self):
        # A chained IndexError carries 'list index out of range'; the
        # StructuralAssemblyError rule must win over the frag_vector rule.
        err = (
            "Generation/Verification failed at UFF_1: StructuralAssemblyError: "
            "Could not assemble a valid 3D complex: every embed attempt failed with "
            "the same structural error (IndexError: list index out of range)\n"
            "Traceback (most recent call last): ..."
        )
        cls, _ = classify_failures.classify(self._report(err))
        self.assertEqual(cls, "structural_assembly")

    def test_generic_pool_exhaustion_still_no_conformers(self):
        err = (
            "Generation/Verification failed at UFF_1: ValueError: "
            "MetalloGen failed to generate any conformers for m-SMILES '[Pt]...'\n"
            "Traceback (most recent call last): ..."
        )
        cls, _ = classify_failures.classify(self._report(err))
        self.assertEqual(cls, "no_conformers")


if __name__ == "__main__":
    unittest.main()
