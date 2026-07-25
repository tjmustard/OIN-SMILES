"""Canonical ligand body (v0.4.5 Lane 1) -- guard the contract.

The graph handed to ``MolToSmiles`` is not canonical even though the serializer is:
perception from 3D distances picks one of several Kekule structures and one of several
resonance forms, so two geometries of the same molecule serialize differently.
``OIN_CANONICAL_BODY`` opts into round-tripping the body through
``MolFromSmiles``/``MolToSmiles`` -- the comparison layer's own fix, promoted upstream.

The contract these guards lock in:

* **default OFF** -> byte-identical output; the new module is not even imported;
* the two forms of a drifting ligand **converge** on one canonical body;
* the ``{n}`` markers still land on the **same donor atoms** -- slot identity is carried
  through the reparse by atom map number, never re-derived;
* every guard **bails to the un-reparsed body for the WHOLE fragment**, never partially;
* the emitted body is a **fixed point** of ``canonical_body``, so it is a canonical
  representative rather than one step along a walk;
* fac and mer stay **distinct** -- folding notation drift must not fold a real isomer.
"""

import os
import sys
import unittest
from pathlib import Path

from rdkit import Chem, RDLogger

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.oin.canonical_body import (  # noqa: E402
    canonical_body,
    canonical_body_emit,
)
from oinsmiles.oin.compare import _SLOT_RE, canonical_roundtrip_key  # noqa: E402
from oinsmiles.oin.inline import OINInlineHandler  # noqa: E402

RDLogger.DisableLog("rdApp.*")

FIX = _ROOT / "tests" / "fixtures"

#: Fast fixtures spanning the awkward cases Lane 1 must not break: eta rings (Ferrocene,
#: TiCat1), a deprotonated X-type aryl carbon (Ir(ppy)3), a chelate (en), a metal-bound
#: secondary amine (POJJOP), phosphines (BINAP/BDPP), carbonyls (FeCO5), a hydride
#: (FeH2(CO)4), and a Kekule/aromatic porphyrin that actually moves (YESKOZ).
_FIXTURES = [
    "CisPlatin.xyz",
    "Cis-PtCl2(en).xyz",
    "Ferrocene.xyz",
    "Ferrocene-halide-face.xyz",
    "fac-Ir(ppy)3.xyz",
    "mer-Ir(ppy)3.xyz",
    "FeCO5.xyz",
    "FeH2(CO)4.xyz",
    "TiCat1.xyz",
    "PdCl2-R-BINAP.xyz",
    "PdCl2-RR-BDPP.xyz",
    "POJJOP.xyz",
    "PtMeNH3ClBr-Cis.xyz",
    "YESKOZ.xyz",
]

#: The two ways a maleic-anhydride-like eta ligand gets perceived -- the exact drift the
#: capstone sweep shows for APABAV_comp_0. Aromatic from the crystal, Kekule from the
#: regenerated structure.
_APABAV_AROMATIC = "O=c1[cH][cH]c(=O)o1"
_APABAV_KEKULE = "O=C1[CH]=[CH]C(=O)O1"


def _lever(name, value):
    """Set an env lever; returns the callable that restores the previous state."""
    prev = os.environ.get(name)

    def restore():
        if prev is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = prev

    if value:
        os.environ[name] = "1"
    else:
        os.environ.pop(name, None)
    return restore


def _flag(value):
    """``OIN_CANONICAL_BODY`` lever; returns its restore callable."""
    return _lever("OIN_CANONICAL_BODY", value)


def _bodies(oin):
    """Slot-stripped ligand bodies of an OIN string, metal fragment excluded."""
    out = []
    for frag in (oin or "").split("."):
        if OINInlineHandler.METAL_REGEX.fullmatch(frag):
            continue
        body = _SLOT_RE.sub("", frag).strip()
        if body:
            out.append(body)
    return out


def _slot_donor_pairs(oin):
    """``{(slot, donor element symbol)}`` for every marker in an OIN string.

    The marker's donor is the atom immediately preceding it, which is what
    ``OINInlineHandler._count_smiles_atoms_before`` resolves at parse time. Comparing this
    multiset between the flag-off and flag-on arms is the marker-safety check: if the
    reparse moved a marker onto a different atom, the element or the slot changes.
    """
    pairs = []
    for frag in (oin or "").split("."):
        if OINInlineHandler.METAL_REGEX.fullmatch(frag):
            continue
        for m in _SLOT_RE.finditer(frag):
            head = frag[: m.start()]
            # Trailing atom token: a bracket atom, or the last organic-subset symbol.
            if head.endswith("]"):
                token = head[head.rindex("[") + 1 :]
                symbol = "".join(ch for ch in token if ch.isalpha())[:2]
            else:
                symbol = head[-2:] if head[-2:] in ("Cl", "Br") else head[-1:]
            pairs.append((m.group(0), symbol.capitalize()))
    return sorted(pairs)


class TestFlagOffIsByteIdentical(unittest.TestCase):
    """Unset lever -> pristine output. The whole opt-in contract rests on this."""

    def setUp(self):
        self.addCleanup(_flag(False))

    def test_goldens_unchanged(self):
        conv = XYZToSMILES()
        self.assertEqual(
            conv.convert(str(FIX / "CisPlatin.xyz")),
            "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
        )
        self.assertEqual(
            conv.convert(str(FIX / "Cis-PtCl2(en).xyz")),
            "[Pt_SPL].[NH2]{0}CC[NH2]{1}.[Cl]{2}.[Cl]{3}",
        )
        self.assertEqual(
            conv.convert(str(FIX / "Ferrocene.xyz")),
            "[Fe_LIN]."
            "[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1."
            "[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
        )

    def test_module_not_imported_when_flag_unset(self):
        """The lazy import is the mechanism that makes 'flag off' cost nothing."""
        import subprocess

        code = (
            "import sys; sys.path.insert(0, %r);"
            "from oinsmiles import XYZToSMILES;"
            "XYZToSMILES().convert(%r);"
            "print('oinsmiles.oin.canonical_body' in sys.modules)"
        ) % (str(_ROOT / "src"), str(FIX / "CisPlatin.xyz"))
        env = {k: v for k, v in os.environ.items() if k != "OIN_CANONICAL_BODY"}
        out = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, env=env
        ).stdout
        self.assertIn("False", out, f"canonical_body must not be imported when off: {out!r}")


class TestCanonicalBodyFunction(unittest.TestCase):
    """``canonical_body`` -- the function Lane 2 calls for its vertex colors."""

    def test_kekule_and_aromatic_converge(self):
        """The APABAV drift: two perceptions of one ligand, one canonical body."""
        a = canonical_body(_APABAV_AROMATIC)
        k = canonical_body(_APABAV_KEKULE)
        self.assertEqual(a, k)
        self.assertFalse(a.startswith("RAW:"), a)

    def test_is_idempotent_on_its_own_output(self):
        for body in (_APABAV_AROMATIC, _APABAV_KEKULE, "NCCN", "[c]1ccccc1-c1ccccn1"):
            once = canonical_body(body)
            self.assertEqual(canonical_body(once), once, body)

    def test_unparseable_body_gets_stable_raw_token(self):
        """A borane cluster / over-valent C#O must still compare by string, not explode."""
        self.assertEqual(canonical_body("C#O"), "RAW:C#O")


class TestCanonicalBodyEmitGuards(unittest.TestCase):
    """Every failure mode must bail, and bail for the WHOLE fragment."""

    @staticmethod
    def _frag(smiles):
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, smiles
        return mol

    def test_returns_fixed_point_of_canonical_body(self):
        for smiles, donors in [
            (_APABAV_KEKULE, [2, 3]),
            (_APABAV_AROMATIC, [1, 2]),
            ("NCCN", [0, 3]),
            ("[O-]C(=O)C", [0]),
            ("O=S(=O)([O-])c1ccccc1", [3]),
            ("CN1C=CN(C)C1", [6]),  # NHC carbene carbon
        ]:
            with self.subTest(smiles=smiles):
                got = canonical_body_emit(self._frag(smiles), donors)
                self.assertIsNotNone(got, smiles)
                body = got[0]
                self.assertEqual(canonical_body(body), body, f"{smiles} -> {body}")

    def test_no_atom_map_residue_in_emitted_body(self):
        got = canonical_body_emit(self._frag("NCCN"), [0, 3])
        self.assertIsNotNone(got)
        self.assertNotIn(":", got[0])

    def test_preexisting_map_number_bails(self):
        """We could not tell our labels from the caller's, so refuse rather than guess."""
        mol = Chem.RWMol(self._frag("NCCN"))
        mol.GetAtomWithIdx(1).SetAtomMapNum(7)
        self.assertIsNone(canonical_body_emit(mol.GetMol(), [0, 3]))

    def test_out_of_range_donor_bails(self):
        self.assertIsNone(canonical_body_emit(self._frag("NCCN"), [0, 99]))

    def test_donor_positions_index_the_donor_element(self):
        """The returned positions must point at atoms of the donor's own element."""
        mol = self._frag("NCCN")
        got = canonical_body_emit(mol, [0, 3])
        self.assertIsNotNone(got)
        body, positions, _reparsed = got
        emitted = Chem.MolFromSmiles(body)
        self.assertIsNotNone(emitted)
        for donor_idx, pos in positions.items():
            self.assertEqual(
                emitted.GetAtomWithIdx(pos).GetAtomicNum(),
                mol.GetAtomWithIdx(donor_idx).GetAtomicNum(),
            )

    def test_oscillating_body_bails_rather_than_emitting_a_non_fixed_point(self):
        """RDKit flips @/@@ on adamantane-cage carbons every parse/write cycle.

        There is no fixed point to reach, so the only correct answer is to keep the
        un-reparsed body. Six of 6062 capstone bodies behave this way. If RDKit ever
        stabilizes this, the assertion below flips to an equality and that is a *win* --
        the point of pinning it is that we notice.
        """
        cage = "CC(C)c1cccc(C(C)C)c1N1C[C@]2(CC1(C)C)[C@H]1C[C@H]3C[C@@H](C1)C[C@@H]2C3"
        once = canonical_body(cage)
        self.assertNotEqual(canonical_body(once), once)
        mol = Chem.MolFromSmiles(cage)
        self.assertIsNotNone(mol)
        donors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N"]
        self.assertIsNone(canonical_body_emit(mol, donors))

    def test_chelate_locked_double_bond_loses_its_ez_marker(self):
        """A C=C locked by a ring that closes through the metal has no free E/Z.

        Both ends of the double bond must sit INSIDE the chelate ring: with N and O as the
        donors, the ring is Fe-N-C=C-O-Fe, so the alkene cannot rotate and its ``/`` is an
        artefact of the SMILES traversal that flips between round-trip directions.
        """
        mol = Chem.MolFromSmiles(r"N/C=C/O")
        self.assertIsNotNone(mol)
        got = canonical_body_emit(mol, [0, 3])
        self.assertIsNotNone(got)
        self.assertNotIn("/", got[0])
        self.assertNotIn("\\", got[0])

    def test_double_bond_only_half_inside_the_metal_ring_is_untouched(self):
        """One end in the ring is not locked -- ``C/C=N/CCO`` still rotates freely."""
        mol = Chem.MolFromSmiles(r"C/C=N/CCO")
        self.assertIsNotNone(mol)
        got = canonical_body_emit(mol, [2, 5])
        self.assertIsNotNone(got)
        self.assertTrue("/" in got[0] or "\\" in got[0], got[0])

    def test_pendant_double_bond_keeps_its_ez_marker(self):
        """A double bond in NO metal ring is genuinely flippable and must be untouched."""
        mol = Chem.MolFromSmiles(r"C/C=C/CCN")
        self.assertIsNotNone(mol)
        got = canonical_body_emit(mol, [5])  # only the N binds; the alkene is pendant
        self.assertIsNotNone(got)
        self.assertTrue("/" in got[0] or "\\" in got[0], got[0])


class TestFlagOnEncoderIntegration(unittest.TestCase):
    """End-to-end: the lever must move the body and nothing else."""

    def setUp(self):
        self.addCleanup(_flag(True))
        self.conv = XYZToSMILES()

    def _pair(self, fixture):
        restore = _flag(False)
        try:
            off = XYZToSMILES().convert(str(FIX / fixture))
        finally:
            restore()
        _flag(True)
        on = XYZToSMILES().convert(str(FIX / fixture))
        return off, on

    def test_markers_stay_on_the_same_donors(self):
        """The safety property: a reparse must never move a {n} onto another atom."""
        for fixture in _FIXTURES:
            with self.subTest(fixture=fixture):
                off, on = self._pair(fixture)
                self.assertEqual(_slot_donor_pairs(off), _slot_donor_pairs(on), f"{off} -> {on}")

    def test_isomer_is_preserved(self):
        """Turning the lever on must not change WHICH isomer the string denotes."""
        for fixture in _FIXTURES:
            with self.subTest(fixture=fixture):
                off, on = self._pair(fixture)
                self.assertEqual(canonical_roundtrip_key(off), canonical_roundtrip_key(on))

    def test_facmer_stays_distinct(self):
        """Over-folding guard: folding notation drift must not fold a real isomer."""
        fac = self.conv.convert(str(FIX / "fac-Ir(ppy)3.xyz"))
        mer = XYZToSMILES().convert(str(FIX / "mer-Ir(ppy)3.xyz"))
        self.assertNotEqual(fac, mer)
        self.assertNotEqual(canonical_roundtrip_key(fac), canonical_roundtrip_key(mer))

    def test_every_fixture_still_encodes(self):
        for fixture in _FIXTURES:
            with self.subTest(fixture=fixture):
                got = self.conv.convert(str(FIX / fixture))
                self.assertTrue(got and "_" in got.split(".")[0], got)
                self.assertTrue(_bodies(got), got)


class TestCanonicalPerception(unittest.TestCase):
    """``OIN_CANONICAL_PERCEPTION`` -- the half of Lane 1 that a reparse cannot do.

    Perception from 3D distances resolves two things by *input atom order*: which of
    several equally-good resonance forms ``lig_checks``/``_select_lig_mol`` keeps (ties go
    to whichever the supplier yielded first), and which valence walk / Kekule matching
    ``AC2BO`` finds first. Renumbering the atoms in the XYZ file therefore hands the
    serializer a different graph -- and re-serializing a different graph cannot repair it,
    which is exactly why the body reparse is measured to be a no-op on this class.
    """

    def setUp(self):
        self.addCleanup(_lever("OIN_CANONICAL_PERCEPTION", False))
        self.addCleanup(_flag(False))

    @staticmethod
    def _renumbered(path, seed):
        """Write a copy of an XYZ with its atom lines permuted. Same molecule, same 3D."""
        import random
        import tempfile

        lines = Path(path).read_text().splitlines()
        n = int(lines[0].split()[0])
        body = lines[2 : 2 + n]
        order = list(range(n))
        random.Random(seed).shuffle(order)
        tmp = Path(tempfile.mkdtemp()) / "renumbered.xyz"
        tmp.write_text("\n".join([lines[0], lines[1]] + [body[i] for i in order]) + "\n")
        return str(tmp)

    def test_atom_permutation_is_renumbering_invariant(self):
        """The permutation must come from the write order, not from CanonicalRankAtoms.

        Measured: over 20 renumberings of ``CC(N)=NC``,
        ``CanonicalRankAtoms(breakTies=True)`` returned a different ranking 18 times -- it
        settles ties between symmetry-equivalent atoms on the input index. ``MolToSmiles``
        returned one string every time. This pins that we build on the invariant one.
        """
        import random

        import numpy as np

        from oinsmiles.utils.xyz2mol_local import _canonical_atom_permutation

        base = Chem.AddHs(Chem.MolFromSmiles("CC(N)=NC"))
        n = base.GetNumAtoms()
        ac = np.zeros((n, n), dtype=int)
        for bond in base.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            ac[i, j] = ac[j, i] = 1
        atoms = [a.GetAtomicNum() for a in base.GetAtoms()]

        def canonical_form(atoms_, ac_):
            perm = _canonical_atom_permutation(ac_, atoms_)
            self.assertIsNotNone(perm)
            idx = np.asarray(perm)
            return tuple(atoms_[i] for i in perm), ac_[np.ix_(idx, idx)].tobytes()

        expected = canonical_form(atoms, ac)
        rng = random.Random(7)
        for _ in range(15):
            p = list(range(n))
            rng.shuffle(p)
            idx = np.asarray(p)
            self.assertEqual(
                canonical_form([atoms[i] for i in p], ac[np.ix_(idx, idx)]),
                expected,
            )

    def test_naxdoi_drifts_under_renumbering_with_the_lever_off(self):
        """The defect this lever exists for. If this ever stops failing, the lever is moot.

        NAXDOI is the smallest in-repo fixture that reproduces it: permuting its atom lines
        changes the emitted string AND the comparison key, because the resonance form the
        tie-break keeps changes with the numbering.
        """
        conv = XYZToSMILES()
        base = conv.convert(str(FIX / "NAXDOI.xyz"))
        drifted = [
            XYZToSMILES().convert(self._renumbered(str(FIX / "NAXDOI.xyz"), seed))
            for seed in (1, 2, 3)
        ]
        self.assertTrue(
            any(d != base for d in drifted),
            "NAXDOI no longer drifts under renumbering -- re-derive the fixture",
        )

    def test_naxdoi_is_renumbering_invariant_with_the_lever_on(self):
        _lever("OIN_CANONICAL_PERCEPTION", True)
        conv = XYZToSMILES()
        base = conv.convert(str(FIX / "NAXDOI.xyz"))
        for seed in (1, 2, 3):
            got = XYZToSMILES().convert(self._renumbered(str(FIX / "NAXDOI.xyz"), seed))
            self.assertEqual(got, base, f"seed {seed}")
            self.assertEqual(canonical_roundtrip_key(got), canonical_roundtrip_key(base))

    def test_lever_off_is_byte_identical(self):
        conv = XYZToSMILES()
        self.assertEqual(
            conv.convert(str(FIX / "CisPlatin.xyz")),
            "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
        )
        self.assertEqual(
            conv.convert(str(FIX / "fac-Ir(ppy)3.xyz")),
            "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1",
        )

    def test_lever_on_preserves_the_isomer_on_the_goldens(self):
        _lever("OIN_CANONICAL_PERCEPTION", True)
        _flag(True)
        for fixture in ("CisPlatin.xyz", "Cis-PtCl2(en).xyz", "Ferrocene.xyz", "FeCO5.xyz"):
            with self.subTest(fixture=fixture):
                restore = _lever("OIN_CANONICAL_PERCEPTION", False)
                off_restore = _flag(False)
                off = XYZToSMILES().convert(str(FIX / fixture))
                restore()
                off_restore()
                _lever("OIN_CANONICAL_PERCEPTION", True)
                _flag(True)
                on = XYZToSMILES().convert(str(FIX / fixture))
                self.assertEqual(canonical_roundtrip_key(off), canonical_roundtrip_key(on))


if __name__ == "__main__":
    unittest.main()
