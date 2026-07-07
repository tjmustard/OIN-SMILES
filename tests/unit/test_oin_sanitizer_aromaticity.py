"""Unit pin for TASK-42 (WS-2, root cause E1): OINSanitizer.generate_robust_smiles
restores aromatic bond type on ring bonds the OIN->XYZ generator left as SINGLE.

The generator de-aromatizes rings for ETKDG (aromatic -> SINGLE) and never
restores them, so re-encoding its bonded mol would serialize an aromatic Cp as
[cH]-[cH]-... . The sanitizer now re-aromatizes ring bonds between aromatic
atoms, while leaving genuine biaryl single bonds (guarded by IsInRing) intact.

The staged SanitizeMol approach the handoff proposed does NOT work here and is
pinned as a negative: a charge-less Cp anion raises KekulizeException on full
sanitize, and SANITIZE_ALL ^ KEKULIZE + SetAromaticity leaves the bond types
SINGLE (dashes survive). See spec/worklog/TASK-42 for the probe evidence.
"""

import unittest

from rdkit import Chem

from oinsmiles.utils.oin_aligner import OINSanitizer


def _broken_cp():
    """A cyclopentadienyl ring in the generator's de-aromatized state:
    5 aromatic-flagged carbons (1 H each) joined by SINGLE bonds."""
    rw = Chem.RWMol()
    idx = []
    for _ in range(5):
        a = Chem.Atom(6)
        a.SetIsAromatic(True)
        a.SetNumExplicitHs(1)
        a.SetNoImplicit(True)
        idx.append(rw.AddAtom(a))
    for i in range(5):
        rw.AddBond(idx[i], idx[(i + 1) % 5], Chem.BondType.SINGLE)
    return rw.GetMol()


def _biaryl_deramatized_rings():
    """Biphenyl with its two benzene rings de-aromatized to SINGLE ring bonds
    but atoms still aromatic-flagged; the inter-ring biaryl bond is a genuine
    (non-ring) SINGLE bond that must be preserved."""
    bph = Chem.MolFromSmiles("c1ccccc1-c1ccccc1")
    rw = Chem.RWMol(bph)
    for b in rw.GetBonds():
        if b.IsInRing() and b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
    return rw.GetMol()


class TestOinSanitizerAromaticity(unittest.TestCase):
    def test_broken_cp_serializes_aromatic_not_single(self):
        """Aromatic-atom + SINGLE-bond Cp -> aromatic SMILES with no explicit
        single bonds between ring atoms (was [cH]-[cH]-[cH]-[cH]-[cH]-)."""
        smi, _ = OINSanitizer.generate_robust_smiles(_broken_cp(), [0])
        # Core regression guard: was [cH]-[cH]-[cH]-[cH]-[cH]-, now dash-free.
        self.assertNotIn("-", smi, f"aromatic ring emitted explicit single bonds: {smi!r}")
        # Parse without sanitize: a neutral cyclopentadienyl is not a valid
        # closed-shell aromatic to RDKit, so default MolFromSmiles would reject
        # it; we only need to confirm the ring bonds are typed AROMATIC.
        m = Chem.MolFromSmiles(smi, sanitize=False)
        self.assertIsNotNone(m, f"output SMILES did not parse: {smi!r}")
        Chem.FastFindRings(m)
        ring_bonds = [b for b in m.GetBonds() if b.IsInRing()]
        self.assertTrue(ring_bonds, "expected a ring")
        self.assertTrue(
            all(b.GetBondType() == Chem.BondType.AROMATIC for b in ring_bonds),
            f"ring bonds not all aromatic: {smi!r}",
        )

    def test_biaryl_single_bond_preserved(self):
        """The IsInRing guard must leave a genuine biaryl single bond intact
        while still re-aromatizing the two rings."""
        smi, _ = OINSanitizer.generate_robust_smiles(_biaryl_deramatized_rings(), [0])
        self.assertIn("-", smi, f"biaryl single bond was wrongly aromatized: {smi!r}")
        m = Chem.MolFromSmiles(smi)
        self.assertIsNotNone(m, f"output SMILES did not parse: {smi!r}")
        # Two separate aromatic six-membered rings joined by one non-ring bond.
        self.assertEqual(m.GetRingInfo().NumRings(), 2, f"expected 2 rings: {smi!r}")
        inter_ring = [b for b in m.GetBonds() if not b.IsInRing()]
        self.assertTrue(
            any(b.GetBondType() == Chem.BondType.SINGLE for b in inter_ring),
            f"biaryl bond not preserved as SINGLE: {smi!r}",
        )


if __name__ == "__main__":
    unittest.main()
