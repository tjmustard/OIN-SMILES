"""Non-binding 0-H chalcogen donors must not be re-protonated on round-trip.

``OINSanitizer.generate_robust_smiles`` locks the hydrogen count of the metal
*binding* atoms it is handed, but historically left every other atom to RDKit's
default implicit-H behaviour. A non-binding heteroatom that carries a genuine
valence deficit -- a croconate/oxo ring O, a nitrito ``-O`` -- has no hydrogen
and no metal bond to fill its valence, yet ``MolToSmiles`` serialized it BARE
(``c(O)``, ``ON=O``): ``SetNoImplicit`` alone does not force a bracket. Any
downstream reader (the MetalloGen adapter's ``MolFromSmiles``) then re-added the
implicit hydrogen, so the regenerated 3D structure gained a phantom H the input
never had -- a round-trip atom-count mismatch (COLWIK croconate 55->58, ACOXEX
nitrito 75->77).

The fix charges such an atom by its deficit (a valence-1 O becomes ``[O-]``):
0 H, closed-shell, unambiguous in the string, and embeddable (a neutral 0-H O
is a radical, which the adapter drops back to bare O and which UFF cannot type).
These tests pin that a bare non-binding chalcogen deficit is closed with a
charge, while a real O-H / carbonyl / ether O is never touched.
"""

import unittest

from rdkit import Chem

from oinsmiles.utils.oin_aligner import OINSanitizer


def _frag_with_bare_deficit_O():
    """Build ``C-O`` where the O mirrors a non-binding deficit donor:
    valence 1, neutral, 0 radicals, 0 H, ``NoImplicit`` set (the exact state a
    croconate/nitrito O reaches ``generate_robust_smiles`` in)."""
    rw = Chem.RWMol()
    c = rw.AddAtom(Chem.Atom(6))
    o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c, o, Chem.BondType.SINGLE)
    c_atom = rw.GetAtomWithIdx(c)
    c_atom.SetNoImplicit(True)
    c_atom.SetNumExplicitHs(3)  # methyl carbon
    o_atom = rw.GetAtomWithIdx(o)
    o_atom.SetNoImplicit(True)
    o_atom.SetNumExplicitHs(0)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol, o


def _reparse_h_on(symbol, smiles):
    """Total H on the first atom of ``symbol`` after a serialize -> parse cycle,
    i.e. what a downstream reader (the generator) actually builds."""
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"unparseable: {smiles}"
    for atom in m.GetAtoms():
        if atom.GetSymbol() == symbol:
            return atom.GetTotalNumHs()
    raise AssertionError(f"no {symbol} in {smiles}")


class TestNonBindingDeficitClosed(unittest.TestCase):
    """A bare non-binding chalcogen deficit is closed with a formal charge."""

    def test_bare_nonbinding_oxygen_gets_charge_not_phantom_h(self):
        mol, o_idx = _frag_with_bare_deficit_O()
        smiles, kmol = OINSanitizer.generate_robust_smiles(mol, [])  # O is non-binding
        # The returned mol carries the closed shell explicitly.
        o_atom = kmol.GetAtomWithIdx(o_idx)
        self.assertEqual(o_atom.GetFormalCharge(), -1)
        self.assertEqual(o_atom.GetTotalNumHs(), 0)
        # And, crucially, a downstream reader rebuilds 0 H (the phantom is gone).
        self.assertIn("[O-]", smiles)
        self.assertEqual(_reparse_h_on("O", smiles), 0)

    def test_nitrito_nonbinding_terminal_oxygen(self):
        # O=N-O with the middle N as the metal binder: the terminal single-bonded
        # O is the non-binding deficit donor (ACOXEX nitrito).
        rw = Chem.RWMol()
        o_term = rw.AddAtom(Chem.Atom(8))
        n = rw.AddAtom(Chem.Atom(7))
        o_dbl = rw.AddAtom(Chem.Atom(8))
        rw.AddBond(o_term, n, Chem.BondType.SINGLE)
        rw.AddBond(n, o_dbl, Chem.BondType.DOUBLE)
        for idx in (o_term, o_dbl):
            a = rw.GetAtomWithIdx(idx)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(0)
        rw.GetAtomWithIdx(n).SetNoImplicit(True)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        smiles, kmol = OINSanitizer.generate_robust_smiles(mol, [n])  # N binds
        self.assertEqual(kmol.GetAtomWithIdx(o_term).GetFormalCharge(), -1)
        self.assertIn("[O-]", smiles)
        # The double-bonded O has full valence: untouched, still 0 H, no charge.
        self.assertEqual(kmol.GetAtomWithIdx(o_dbl).GetFormalCharge(), 0)


class TestRealHydrogensPreserved(unittest.TestCase):
    """Deficit-free / hydrogen-bearing donors are never charged."""

    def test_real_hydroxyl_keeps_its_hydrogen(self):
        # Methanol: the O has an implicit hydrogen (full valence) and must stay
        # a neutral O-H, never become [O-].
        mol = Chem.MolFromSmiles("CO")
        smiles, kmol = OINSanitizer.generate_robust_smiles(mol, [])
        o_atom = next(a for a in kmol.GetAtoms() if a.GetSymbol() == "O")
        self.assertEqual(o_atom.GetFormalCharge(), 0)
        self.assertNotIn("[O-]", smiles)
        self.assertGreaterEqual(_reparse_h_on("O", smiles), 1)

    def test_carbonyl_and_ether_oxygen_untouched(self):
        # Acetone (C=O, valence 2) and dimethyl ether (C-O-C, valence 2): no
        # deficit, so no charge is added.
        for smi in ("CC(C)=O", "COC"):
            mol = Chem.MolFromSmiles(smi)
            smiles, kmol = OINSanitizer.generate_robust_smiles(mol, [])
            self.assertNotIn("[O-]", smiles, f"{smi} wrongly charged")
            for a in kmol.GetAtoms():
                if a.GetSymbol() == "O":
                    self.assertEqual(a.GetFormalCharge(), 0)

    def test_nitrogen_is_not_charged(self):
        # A bare non-binding N (amido/nitride/ammine ambiguity) is out of scope:
        # it must NOT be charged here (owned by oin/inline.py notation).
        mol, _ = _frag_with_bare_deficit_O()
        # swap the O for an N deficit: C#N-like terminal handled elsewhere
        rw = Chem.RWMol()
        c = rw.AddAtom(Chem.Atom(6))
        natom = rw.AddAtom(Chem.Atom(7))
        rw.AddBond(c, natom, Chem.BondType.DOUBLE)
        a = rw.GetAtomWithIdx(natom)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        rw.GetAtomWithIdx(c).SetNoImplicit(True)
        rw.GetAtomWithIdx(c).SetNumExplicitHs(2)
        m = rw.GetMol()
        m.UpdatePropertyCache(strict=False)
        _, kmol = OINSanitizer.generate_robust_smiles(m, [])
        for atom in kmol.GetAtoms():
            if atom.GetSymbol() == "N":
                self.assertEqual(atom.GetFormalCharge(), 0)


if __name__ == "__main__":
    unittest.main()
