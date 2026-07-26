"""Regression fixtures for the v0.4.5 `atom_count` hard-fail class.

That class -- 74 capstone molecules, the largest non-timeout `hard_fail` group -- is
**entirely hydrogen**: over the 27 molecules with a stored generated structure, 27/27
differed in H and nothing else. Two independent causes, one fixture set each. See
``docs/ATOM_COUNT_v0.4.5.md``.

The molecules and OIN strings below are real capstone rows, not constructions.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem, RDLogger

from oinsmiles.generation.metallogen_adapter import _prepare_ligand_fragments
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.oin.hydrogen import h_faithful_smiles

RDLogger.DisableLog("rdApp.*")


def _parse(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
    )
    return mol


def _atom_total(smiles: str) -> int:
    """Atoms a SMILES implies, hydrogens included."""
    mol = _parse(smiles)
    assert mol is not None, f"unparseable: {smiles}"
    return sum(1 + (0 if a.GetAtomicNum() == 1 else a.GetTotalNumHs()) for a in mol.GetAtoms())


def _adapter_atom_total(oin: str) -> int:
    """Atoms the generator would be asked to build from `oin`."""
    parsed = OINParser().parse(oin)
    metal_frag, specs, _geo = _prepare_ligand_fragments(parsed)
    import re

    total = _atom_total(re.sub(r"_[A-Z0-9]+", "", metal_frag))
    for smi, _winding in specs:
        total += _atom_total(smi)
    return total


class TestKekulizeRescueKeepsHydrogen(unittest.TestCase):
    """Cause B: the Cp kekulization rescue must not destroy an innocent aromatic C-H.

    ``_prepare_ligand_fragments`` charges ``ring[0]`` of every 5-membered all-aromatic
    ring when a fragment will not kekulize. A -1 on a BARE aromatic carbon flips its
    implicit-H count from 1 to 0, so a thiophene riding along on a ligand whose carbene
    donor is the real problem silently loses a hydrogen. This accounted for 13 of the
    13 auditable atom-count LOSS rows.
    """

    #: (molecule, OIN string as the capstone sweep recorded it, input XYZ atom count)
    CASES = [
        (
            "QOBFOF_comp_0",
            "[Au_SPL].C[C@]1(O)C{0}=C(c2ccsc2)n2ccccc21.[Cl]{1}.[Cl]{2}.[Cl]{3}",
            31,
        ),
        (
            "AJODEI_comp_0",
            "[Zn_TET]."
            "CCCCCc1cc(C=C2N{0}=C(c3ccc[nH]3)C=C2OC)n{1}c1C."
            "CCCCCc1cc(C=C2N{2}=C(c3ccc[nH]3)C=C2OC)n{3}c1C",
            97,
        ),
        (
            "MUXKOH_comp_0",
            "[Pt_SPL]."
            "Cc1ccc(-c2cc(-c3cccc(-c4cc(-c5ccc(C)cc5)nn{0}4)n{1}3)n{2}n2)cc1."
            "c1ccc(P{3}(c2ccccc2)c2ccccc2)cc1",
            84,
        ),
        (
            "DUJPIJ_comp_0",
            "[Ti_OCT].CN{0}C.c1ccc2n{3}c(-n3cccn{1}3)cc2c1.c1ccc2n{2}c(-n3cccn{4}3)cc2c1.CN{5}C",
            63,
        ),
        (
            "EGUDAL_comp_0",
            "[Ta_OCT].c1ccc(C(c2ccccc2)(c2cccn{3}2)c2cccn{0}2)cc1.[Cl]{1}."
            "c1ccc(C(c2ccccc2)(c2cccn{2}2)c2cccn{4}2)cc1.[Cl]{5}",
            81,
        ),
    ]

    def test_adapter_preserves_input_atom_count(self):
        for name, oin, expected in self.CASES:
            with self.subTest(molecule=name):
                self.assertEqual(
                    _adapter_atom_total(oin),
                    expected,
                    f"{name}: the generator would be asked to build a different number "
                    "of atoms than the input has -- the kekulization rescue is eating a "
                    "hydrogen again",
                )

    def test_cp_ring_still_gets_its_charge_and_keeps_its_hydrogens(self):
        """The rescue must keep working for the case it exists for.

        A Cp arrives as ``[cH]`` brackets with ``NoImplicit`` already set, so freezing
        the H count before charging is a no-op there: the -1 still lands and all five
        C-H survive, exactly as before. ``UQUXAG_comp_0``, a real capstone row that
        carries an unsubstituted Cp.
        """
        oin = (
            "[Fe_TET].[CH2]{0}C[SH]=C1[CH][CH][CH][CH][CH]1."
            "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.C{2}#O.C{3}#O"
        )
        parsed = OINParser().parse(oin)
        _metal, specs, _geo = _prepare_ligand_fragments(parsed)
        cp = [s for s, _w in specs if s.count("c") >= 5 and "Si" not in s]
        self.assertTrue(cp, f"no Cp fragment among {[s for s, _ in specs]}")
        for smi in cp:
            mol = _parse(smi)
            self.assertIsNotNone(mol)
            self.assertTrue(
                any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()),
                f"Cp lost its anionic charge: {smi}",
            )
            # 5 carbons + 5 hydrogens: the charge must not have cost one.
            self.assertEqual(_atom_total(smi), 10, smi)


class TestHFaithfulSmiles(unittest.TestCase):
    """Cause A: a serialized fragment must re-read with the H count it was written with.

    RDKit's writer emits a bare organic-subset symbol whenever it judges brackets
    unnecessary, and a bare symbol re-reads as "fill to the next allowed valence with
    hydrogen" -- ``SetNoImplicit(True)`` does not force a bracket.
    """

    #: fragments whose bare serialization gains hydrogen on re-reading.
    #: (fragment SMILES, atoms it should imply)
    DRIFTING = [
        # CIDDAU_comp_0's thiophene sulfur: perceived valence 3, strictly between
        # sulfur's allowed 2 and 4, so a bare `S` re-reads as `[SH]`.
        ("N#C[C@H]1C=CC=S1", 3),
        # INENOF_comp_0: an NHC carbene donor and a benzylic carbon, both written
        # `[C]` by the encoder and both de-bracketed by a plain re-serialization.
        ("CCCCN1[C]N([C]c2ccccc2)C=C1", 0),
    ]

    def setUp(self):
        self._saved = os.environ.get("OIN_H_FAITHFUL")

    def tearDown(self):
        if self._saved is None:
            os.environ.pop("OIN_H_FAITHFUL", None)
        else:
            os.environ["OIN_H_FAITHFUL"] = self._saved

    @staticmethod
    def _intent_mol(smiles: str):
        """The mol as the encoder holds it: every H count frozen, nothing implicit."""
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        rw = Chem.RWMol(mol)
        for atom in rw.GetAtoms():
            atom.SetNumExplicitHs(atom.GetNumExplicitHs())
            atom.SetNoImplicit(True)
        out = rw.GetMol()
        out.UpdatePropertyCache(strict=False)
        return out

    @staticmethod
    def _total_h(mol) -> int:
        return sum(1 if a.GetAtomicNum() == 1 else a.GetTotalNumHs() for a in mol.GetAtoms())

    def test_lever_off_is_byte_identical_to_moltosmiles(self):
        os.environ.pop("OIN_H_FAITHFUL", None)
        for smiles, _n in self.DRIFTING:
            with self.subTest(fragment=smiles):
                mol = self._intent_mol(smiles)
                self.assertEqual(
                    h_faithful_smiles(mol, isomericSmiles=True, canonical=True),
                    Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
                )

    def test_lever_on_makes_the_string_reread_with_the_same_hydrogen(self):
        os.environ["OIN_H_FAITHFUL"] = "1"
        for smiles, _n in self.DRIFTING:
            with self.subTest(fragment=smiles):
                mol = self._intent_mol(smiles)
                intent = self._total_h(mol)
                plain = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                fixed = h_faithful_smiles(mol, isomericSmiles=True, canonical=True)
                # The premise: the plain string really does drift. If RDKit ever stops
                # doing this, this assertion is how we find out.
                self.assertNotEqual(
                    self._total_h(_parse(plain)),
                    intent,
                    f"{smiles}: plain MolToSmiles no longer drifts -- re-check whether "
                    "this lever is still needed",
                )
                self.assertEqual(
                    self._total_h(_parse(fixed)),
                    intent,
                    f"{smiles}: repaired string still re-reads with the wrong H count ({fixed})",
                )

    def test_never_mutates_the_caller_mol(self):
        os.environ["OIN_H_FAITHFUL"] = "1"
        mol = self._intent_mol(self.DRIFTING[0][0])
        before = [
            (a.GetNumRadicalElectrons(), a.GetNumExplicitHs(), a.GetNoImplicit())
            for a in mol.GetAtoms()
        ]
        h_faithful_smiles(mol, isomericSmiles=True, canonical=True)
        after = [
            (a.GetNumRadicalElectrons(), a.GetNumExplicitHs(), a.GetNoImplicit())
            for a in mol.GetAtoms()
        ]
        self.assertEqual(before, after)

    def test_never_raises_on_a_mol_with_no_property_cache(self):
        """``oin/inline.py`` passes a ``sanitize=False`` parse with no property cache.

        ``GetTotalNumHs()`` on such an atom raises a RDKit pre-condition violation
        rather than returning a wrong answer -- and inline.py wraps its caller in a
        bare ``except Exception`` that silently reroutes to a different slot-tagging
        strategy, so an exception here would not surface as an error at all.
        """
        os.environ["OIN_H_FAITHFUL"] = "1"
        raw = Chem.MolFromSmiles("CCCCN1[C]N([C]c2ccccc2)C=C1", sanitize=False)
        out = h_faithful_smiles(raw, canonical=False)
        self.assertIsInstance(out, str)
        self.assertTrue(out)


if __name__ == "__main__":
    unittest.main()
