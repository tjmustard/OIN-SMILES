"""Unit tests for aromatic perception in the XYZ->OIN re-encode path.

``get_oin_string`` is reached two ways: from ``XYZToSMILES.convert``, where
``get_tmc_mol``'s final ``SanitizeMol`` leaves ring bonds AROMATIC-typed, and from
the round-trip harness, which feeds it MetalloGen's contract mol directly. The
contract mol is kekulized in place (``Chem.Kekulize`` keeps the aromatic flags on),
so its ring bonds are SINGLE/DOUBLE-typed with ``GetIsAromatic() == True``.

The fragment rebuild inside ``get_oin_string`` used ``RWMol.AddBond(u, v, type)``,
which copies the bond TYPE but creates the bond with ``IsAromatic == False`` --
while ``AddAtom`` *does* copy the atom's aromatic flag. ``OINSanitizer`` then
upgraded only the SINGLE ring bonds back to AROMATIC, leaving the doubles, and the
fragment serialized as ``Cc1c=cc(C)=cc=1``: lowercase aromatic atoms carrying
explicit double bonds. That string does not re-parse, so ``oin/compare.py`` falls
back to a ``RAW:`` token and the round-trip key can never match (the
``garbled_aromatic`` bucket).

The encoder must therefore be *representation-neutral*: the same mol, whether its
rings are given aromatic-typed or kekulized, must produce the same OIN string. The
Kekule bond orders also perturb ``CanonicalRankAtoms``, which drives eta
heading/winding selection, so a flag-only fix is not enough -- an aromatic-flagged
bond is normalized to the AROMATIC bond type.
"""

import os
import re
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem, RDLogger

from oinsmiles.core.translator import XYZToSMILES
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.utils.perception_tmc import _repair_mixed_aromaticity, get_oin_string, get_tmc_mol

RDLogger.DisableLog("rdApp.*")

_FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "../fixtures"))

# Aromatic ligands of every flavour the encoder handles: sigma-bound carbanion +
# N donor (ppy), eta5 rings whose winding is canonical-rank sensitive (Cp), and a
# non-aromatic control that must not move.
_AROMATIC_FIXTURES = [
    "fac-Ir(ppy)3",
    "mer-Ir(ppy)3",
    "Ferrocene",
    "TiCat1",
    "TiCat3",
    "TiCat4",
]
_CONTROL_FIXTURES = ["Cis-PtCl2(en)", "TransPlatin"]

# An aromatic atom double-bonded to another aromatic atom ("c=c", "c=1"). Legal
# aromatic SMILES never writes this; a quinoid exocyclic "c(=O)" does not match.
_GARBLED = re.compile(r"[bcnops]\d*=[bcnops]|[bcnops]=\d")


def _encode_both_ways(name):
    """Return (aromatic-typed encode, kekulized encode) of the same fixture."""
    tmc_mol, coords = get_tmc_mol(os.path.join(_FIXTURES, f"{name}.xyz"), 0)
    aromatic = get_oin_string(tmc_mol, coords)
    # Exactly what build_contract_mol produces: Kekule bond types, flags intact.
    kekulized = Chem.RWMol(tmc_mol)
    Chem.Kekulize(kekulized)
    return aromatic, get_oin_string(kekulized.GetMol(), coords)


class TestKekuleInvariance(unittest.TestCase):
    """The encoder must not care how its input mol represents an aromatic ring."""

    def test_aromatic_fixtures_encode_identically(self):
        for name in _AROMATIC_FIXTURES:
            with self.subTest(fixture=name):
                aromatic, kekulized = _encode_both_ways(name)
                self.assertEqual(aromatic, kekulized)

    def test_control_fixtures_unaffected(self):
        for name in _CONTROL_FIXTURES:
            with self.subTest(fixture=name):
                aromatic, kekulized = _encode_both_ways(name)
                self.assertEqual(aromatic, kekulized)

    def test_canonical_roundtrip_keys_match(self):
        # The property the harness actually checks. Redundant with byte equality
        # today, but it is the contract that must not regress.
        for name in _AROMATIC_FIXTURES:
            with self.subTest(fixture=name):
                aromatic, kekulized = _encode_both_ways(name)
                self.assertEqual(
                    canonical_roundtrip_key(aromatic), canonical_roundtrip_key(kekulized)
                )


class TestNoGarbledAromaticNotation(unittest.TestCase):
    """A kekulized input must never emit `c=c` / `c=1` mixed notation."""

    def test_kekulized_input_emits_clean_aromatic_smiles(self):
        for name in _AROMATIC_FIXTURES:
            with self.subTest(fixture=name):
                _, kekulized = _encode_both_ways(name)
                match = _GARBLED.search(kekulized)
                self.assertIsNone(
                    match,
                    f"{name} re-encoded with mixed aromatic/double notation: {kekulized}",
                )

    def test_neutral_arene_fragment_reparses(self):
        # Ir(ppy)3's phenyl is a carbanion donor, so its slot-stripped fragment
        # legitimately fails to re-parse. Ferrocene's Cp is the same. Use the
        # pyridine ring, which is a neutral N donor: after stripping its slot the
        # fragment must be a valid molecule, which is what compare.py needs to
        # avoid a RAW: token.
        _, kekulized = _encode_both_ways("fac-Ir(ppy)3")
        for fragment in kekulized.split("."):
            clean = re.sub(r"\{\d+[<>^]?\}", "", fragment)
            if "n" not in clean:
                continue
            # Sanitization is what compare.py runs; it must at least kekulize.
            mol = Chem.MolFromSmiles(clean, sanitize=False)
            self.assertIsNotNone(mol, f"unparseable fragment: {clean}")
            mol.UpdatePropertyCache(strict=False)
            Chem.SetAromaticity(mol, Chem.AROMATICITY_DEFAULT)
            Chem.Kekulize(Chem.RWMol(mol), clearAromaticFlags=True)


class TestPorphyrinMacrocycle(unittest.TestCase):
    """Pins the macrocycle_perception diagnosis (see docs/KNOWN_LIMITATIONS.md).

    A metalloporphyrin's aromatic macrocycle is perceived only because get_tmc_mol
    sanitizes with the metal attached AND two pyrrole nitrogens carry a -1 charge.
    OIN carries no charges, so the generator builds all four nitrogens neutral and
    the contract mol cannot sanitize -- which is why the re-encode localizes. The
    fix belongs in the donor-charge layer, not in aromatic perception.
    """

    def test_encoder_perceives_the_aromatic_macrocycle(self):
        oin = XYZToSMILES().convert(os.path.join(_FIXTURES, "BEGLUU.xyz"))
        self.assertIn("c1c2n{0}c(", oin)

    def test_forward_encode_is_stable(self):
        path = os.path.join(_FIXTURES, "BEGLUU.xyz")
        first = XYZToSMILES().convert(path)
        second = XYZToSMILES().convert(path)
        self.assertEqual(first, second)

    def test_two_pyrrole_nitrogens_are_anionic(self):
        # The charge the OIN string cannot carry. If this ever becomes 0, the
        # macrocycle will stop being perceived aromatic and smiles_1 will localize.
        tmc_mol, _ = get_tmc_mol(os.path.join(_FIXTURES, "BEGLUU.xyz"), 0)
        metal = next(a for a in tmc_mol.GetAtoms() if a.GetSymbol() == "Ni")
        charges = sorted(n.GetFormalCharge() for n in metal.GetNeighbors())
        self.assertEqual(charges, [-1, -1, 0, 0])


class TestRepairMixedAromaticity(unittest.TestCase):
    """`_repair_mixed_aromaticity` re-perceives only when flags and orders disagree."""

    @staticmethod
    def _kekule_without_bond_flags(smiles):
        """Aromatic atoms, Kekule bond types, no aromatic bond flags at all.

        What a producer that rebuilds ring bonds from scratch hands us (the bond
        flag is not recoverable from the bond type alone).
        """
        rw = Chem.RWMol(Chem.MolFromSmiles(smiles))
        Chem.Kekulize(rw, clearAromaticFlags=False)
        for bond in rw.GetBonds():
            bond.SetIsAromatic(False)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        return mol

    def test_repairs_flagless_kekule_ring(self):
        broken = self._kekule_without_bond_flags("Cc1ccc(C)cc1")
        repaired = _repair_mixed_aromaticity(broken)
        ring_bonds = [b for b in repaired.GetBonds() if b.IsInRing()]
        self.assertTrue(ring_bonds)
        self.assertTrue(all(b.GetIsAromatic() for b in ring_bonds))
        self.assertEqual(Chem.MolToSmiles(repaired), "Cc1ccc(C)cc1")

    def test_noop_on_dearomatized_cp(self):
        # The OIN->XYZ generator hands ETKDG an all-SINGLE Cp with aromatic atoms.
        # OINSanitizer restores those rings; the repair must not touch them.
        rw = Chem.RWMol(Chem.MolFromSmiles("c1cc[cH-]c1"))
        for atom in rw.GetAtoms():
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
        for bond in rw.GetBonds():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
        cp = rw.GetMol()
        cp.UpdatePropertyCache(strict=False)
        self.assertIs(_repair_mixed_aromaticity(cp), cp)

    def test_noop_on_aromatic_typed_ring(self):
        benzene = Chem.MolFromSmiles("c1ccccc1")
        self.assertIs(_repair_mixed_aromaticity(benzene), benzene)


if __name__ == "__main__":
    unittest.main()
