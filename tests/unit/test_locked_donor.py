"""Metal-locked donor chirality (Y1 P3) -- guard the contract.

``POJJOP``'s only stereocentre is a Pd-bound secondary amine, and it used to encode
byte-identically to its mirror: total encoder blindness. ``OIN_EMIT_LOCKED_DONOR`` opts
into a descriptor that separates them. The contract locked in here:

* **default OFF** -> byte-identical, so nothing in the corpus moves until it is promoted;
* **lever ON** -> the two enantiomers' raw strings diverge, and the difference is confined
  to chiral tags (the mirror is the base string with every tag swapped);
* the descriptor is a genuine chirality descriptor, not merely a distinguisher: invariant
  under input atom renumbering, invariant under proper rotation, inverting under
  reflection;
* a metal-bound **ammine** (M,H,H,H) and a **primary** amine (M,H,H,R) emit **nothing**,
  because they are not stereocentres -- the over-sensitivity failure the axial lane had to
  fix once already.

The mirror assertions deliberately compare **whole strings under an ``@``<->``@@`` swap**
rather than counting tags. Counting is blind to a symmetric swap: three ``@@`` plus three
``@`` mirrors to three ``@`` plus three ``@@``, identical counts. And a descriptor made
stable by being made *constant* would pass every check written against the easy fixture --
which is exactly what the Y2 axial wave shipped.
"""

import os
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
for _p in (str(_ROOT), str(_ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from rdkit import Chem  # noqa: E402

from oinsmiles import XYZToSMILES  # noqa: E402
from oinsmiles.core.chirality import CIPAssigner  # noqa: E402
from oinsmiles.oin.compare import (  # noqa: E402
    normalize_oin_for_comparison,
    winding_canonical_key,
)
from oinsmiles.oin.locked_donor import ENV_LEVER, plan_locked_donors  # noqa: E402
from oinsmiles.utils.xyz2mol import get_tmc_mol  # noqa: E402

FIX = _ROOT / "tests" / "fixtures"

#: Sole stereocentre is a Pd-bound 2 degree amine -- the motivating P3 case.
POJJOP = str(FIX / "POJJOP.xyz")
#: Three Cu-bound 2 degree amines on one cyclohexane backbone -- the multi-centre case,
#: and the reason the mirror check compares whole strings rather than one tag.
RIFGUJ = str(FIX / "RIFGUJ_comp_2.xyz")
#: Over-sensitivity controls.
CISPLATIN = str(FIX / "CisPlatin.xyz")  # Pt(NH3)2Cl2 -- ammine, M,H,H,H
PT_AMMINE = str(FIX / "PtMeNH3ClBr-Cis.xyz")  # a second ammine
PTCL2_EN = str(FIX / "Cis-PtCl2(en).xyz")  # ethylenediamine -- primary amine, M,H,H,C
BINAP = str(FIX / "PdCl2-R-BINAP.xyz")  # P donors, no locked N


def _swap_tags(s: str) -> str:
    """Exchange every ``@`` and ``@@`` in *s*."""
    return s.replace("@@", "\x00").replace("@", "@@").replace("\x00", "@")


def _swap_n_tags(s: str) -> str:
    """Exchange ``[N@...]`` <-> ``[N@@...]`` and leave every other tag alone.

    The strict "mirror == swap(every tag)" form is wrong in general -- Lane 8 recorded the
    same thing for an eta complex whose mirror differed only in the winding character. Here
    the counter-example is a 1,3,5-trisubstituted cyclohexane backbone (``RIFGUJ``): its
    ring-carbon tags encode a *relative* (all-cis) arrangement and are measured NOT to
    invert under reflection, with the lever off as well as on. So the mirror assertion is
    "every locked-nitrogen descriptor inverted and nothing else in the string moved", which
    is still a whole-string comparison and still immune to the count-blindness that made an
    earlier hand-run of this check report a false pass.
    """
    out = s.replace("[N@@", "\x00").replace("[N@", "[N@@")
    return out.replace("\x00", "[N@")


def _n_descriptors(s: str) -> list[str]:
    """The locked-nitrogen descriptors, in the order they appear in *s*."""
    return re.findall(r"\[N@@?H?\]", s)


def _read_atoms(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    rows = [ln.split() for ln in lines[2 : 2 + n]]
    return lines[1], [(r[0], np.array([float(r[1]), float(r[2]), float(r[3])])) for r in rows]


def _write_atoms(dst: Path, head, atoms) -> str:
    body = "".join(f"{s} {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}\n" for s, c in atoms)
    dst.write_text(f"{len(atoms)}\n{head}\n{body}")
    return str(dst)


class _LeverBase(unittest.TestCase):
    """Set/restore ``OIN_EMIT_LOCKED_DONOR`` around each test.

    ⚠ Pins ``OIN_CANONICAL_BODY`` **OFF**, because the two levers are still INCOMPATIBLE and
    ``OIN_CANONICAL_BODY`` is default-ON since v0.4.5. ``canonical_body_emit`` reparses the
    ligand body, and sanitizing a *metal-free* fragment runs
    ``AssignStereochemistry(cleanIt=True)``, which strips the chiral tag off a 2-degree amine as
    a freely inverting nitrogen -- the exact RDKit behaviour this descriptor exists to work
    around. With both on, the tag is stamped upstream and then discarded, so every assertion
    below would see a bare ``N``.

    **So P3 is not usable in the shipped default configuration.** Stated here rather than hidden,
    because the pin makes this suite green against a configuration nobody ships.

    ⚠ THE OBVIOUS FIX WAS TRIED IN v0.4.6 AND IS MEASURABLY WRONG. Copying the chiral tag onto
    the reparsed donor (the correspondence is available -- ``_reparse_once``'s Guard 2 already
    proves same element, same heavy degree) does make P3 emit under ``OIN_CANONICAL_BODY``, and
    POJJOP passes. But setting a tag after the sanitize introduces a stereocentre the canonical
    ranker did not account for, which changes the canonical WRITE ORDER -- and ``@``/``@@`` is a
    parity relative to that order. On ``RIFGUJ_comp_2`` the three ring-CARBON tags then flip
    between a structure and its mirror, and the geometry says they must not:
    ``AssignStereochemistryFrom3D`` + ``rdCIPLabeler`` label those carbons lowercase ``s``
    (pseudo-asymmetric, a RELATIVE all-cis descriptor) and read ``s`` identically for the
    structure and its reflection.

    ``TestMultiCentreDescriptor::test_flips_under_reflection`` and
    ``TestLeverOnDivergesOnEnantiomers::test_three_locked_amines_all_invert_together`` are what
    caught it. Single-centre POJJOP could not -- the Y2 lesson, intact.

    A correct fix must preserve the tag WITHOUT perturbing the ranking: keep the donor bracketed
    through the sanitize, or re-derive parity from the parent geometry once the write order is
    fixed. See ``oin/canonical_body.py::_reparse_once`` and ``oin/levers.py::_HELD_OFF``.
    """

    def setUp(self):
        self._prev = os.environ.get(ENV_LEVER)
        self.addCleanup(self._restore)
        body = mock.patch.dict(os.environ, {"OIN_CANONICAL_BODY": "0"})
        body.start()
        self.addCleanup(body.stop)

    def _restore(self):
        if self._prev is None:
            os.environ.pop(ENV_LEVER, None)
        else:
            os.environ[ENV_LEVER] = self._prev

    def _set(self, on: bool):
        # "0" rather than unset: OIN_EMIT_LOCKED_DONOR is still default-OFF so both spellings
        # agree today, but being explicit is what stops this from silently inverting if the
        # lever is ever promoted -- the trap that hit five other test modules in v0.4.5.
        os.environ[ENV_LEVER] = "1" if on else "0"

    def _encode(self, path: str) -> str:
        return XYZToSMILES().convert(path)


class TestDefaultOff(_LeverBase):
    """With the lever unset the encoder must behave exactly as it did before."""

    def test_pojjop_blind_by_default(self):
        self._set(False)
        with tempfile.TemporaryDirectory() as d:
            head, atoms = _read_atoms(POJJOP)
            mir = _write_atoms(
                Path(d) / "mirror.xyz", head, [(s, np.array([c[0], c[1], -c[2]])) for s, c in atoms]
            )
            base = self._encode(POJJOP)
            mirror = self._encode(mir)
        self.assertEqual(base, mirror, "default OFF must stay byte-identical (no regression)")
        self.assertIn("[NH]{0}", base, "the bare, configuration-free amine is the OFF behaviour")

    def test_no_locked_tag_appears_when_off(self):
        self._set(False)
        self.assertNotIn("[N@", self._encode(RIFGUJ))


class TestLeverOnDivergesOnEnantiomers(_LeverBase):
    """The point of the lane: the two enantiomers stop colliding."""

    def _base_and_mirror(self, path):
        with tempfile.TemporaryDirectory() as d:
            head, atoms = _read_atoms(path)
            mir = _write_atoms(
                Path(d) / "mirror.xyz", head, [(s, np.array([c[0], c[1], -c[2]])) for s, c in atoms]
            )
            return self._encode(path), self._encode(mir)

    def test_pojjop_enantiomers_diverge(self):
        self._set(True)
        base, mirror = self._base_and_mirror(POJJOP)
        self.assertNotEqual(base, mirror, "P3: the metal-locked amine enantiomers must differ")
        self.assertIn("[N@", base)
        # The difference must be confined to configuration. Comparing under the tag swap
        # (rather than counting tags, which a symmetric swap defeats) also proves no
        # constitution change was smuggled in.
        self.assertEqual(mirror, _swap_tags(base))

    def test_three_locked_amines_all_invert_together(self):
        """RIFGUJ carries three Cu-bound 2 degree amines; every one must invert.

        This is the case a tag *count* comparison cannot police, and the case where a
        descriptor that had been made constant would look fine on POJJOP.
        """
        self._set(True)
        base, mirror = self._base_and_mirror(RIFGUJ)
        self.assertEqual(len(_n_descriptors(base)), 3)
        self.assertNotEqual(base, mirror)
        # Whole-string, not counts: every locked N inverted and nothing else moved.
        self.assertEqual(mirror, _swap_n_tags(base))
        self.assertEqual(set(_n_descriptors(base)) & set(_n_descriptors(mirror)), set())

    def test_round_trip_key_still_folds_the_descriptor(self):
        """The batch harness's acceptance predicate must be blind to this descriptor.

        The key re-parses each fragment through RDKit, which drops trivalent-nitrogen
        chirality, so POJJOP and its mirror still share a key with the lever ON. Two
        consequences, and both are load-bearing: promoting this lever cannot move
        ``facmer_divergent`` or any other harness count, and equally the harness cannot
        *confirm* the descriptor -- only the raw string can.
        """
        self._set(True)
        base, mirror = self._base_and_mirror(POJJOP)
        self.assertNotEqual(base, mirror)
        kb = winding_canonical_key(normalize_oin_for_comparison(base))
        km = winding_canonical_key(normalize_oin_for_comparison(mirror))
        self.assertEqual(kb, km, "the round-trip key must fold the locked-donor descriptor")


class TestDescriptorIsCanonical(_LeverBase):
    """Three-property test on POJJOP: renumbering-, rotation-invariant, reflection-flipping.

    POJJOP carries the descriptor and nothing else stereogenic, and its whole string is
    measured byte-stable under renumbering with the lever off, so the strictest possible
    form applies here: full-string equality. The multi-centre case is
    :class:`TestMultiCentreDescriptor`, which has to be scoped more narrowly for reasons
    recorded there.
    """

    def _encode_variant(self, atoms, head):
        with tempfile.TemporaryDirectory() as d:
            return self._encode(_write_atoms(Path(d) / "probe.xyz", head, atoms))

    def test_invariant_under_atom_renumbering(self):
        self._set(True)
        head, atoms = _read_atoms(POJJOP)
        base = self._encode_variant(atoms, head)
        self.assertIn("[N@", base, "fixture must emit a locked-donor descriptor")
        rng = np.random.default_rng(42)
        for _ in range(3):
            perm = [int(i) for i in rng.permutation(len(atoms))]
            got = self._encode_variant([atoms[i] for i in perm], head)
            self.assertEqual(got, base, "descriptor drifted under input renumbering")

    def test_invariant_under_proper_rotation(self):
        self._set(True)
        head, atoms = _read_atoms(POJJOP)
        base = self._encode_variant(atoms, head)
        rng = np.random.default_rng(7)
        for _ in range(3):
            q, r = np.linalg.qr(rng.normal(size=(3, 3)))
            q = q @ np.diag(np.sign(np.diag(r)))
            if np.linalg.det(q) < 0:
                q[:, 0] *= -1  # force a PROPER rotation: no mirroring
            got = self._encode_variant([(s, q @ c) for s, c in atoms], head)
            self.assertEqual(got, base, "descriptor drifted under proper rotation")

    def test_flips_under_reflection(self):
        self._set(True)
        head, atoms = _read_atoms(POJJOP)
        base = self._encode_variant(atoms, head)
        mirrored = self._encode_variant(
            [(s, np.array([c[0], c[1], -c[2]])) for s, c in atoms], head
        )
        self.assertNotEqual(mirrored, base, "a reflection MUST change the descriptor")
        self.assertEqual(mirrored, _swap_tags(base))


class TestMultiCentreDescriptor(_LeverBase):
    """The three-property test on the three-centre fixture, scoped to what this lane owns.

    ``RIFGUJ``'s whole string is **not** byte-stable under renumbering, and that has
    nothing to do with this lane: measured with the lever **off**, it drifts in 3 of 3
    renumberings -- ring-carbon ``[C@H]``/``[C@@H]`` tags flip, one loses its tag
    entirely, and the ``{2}``/``{3}`` slot numbers swap. That is the pre-existing 13%
    stereo-flip class (Lane 8's ``OIN_STABLE_STEREO``) plus slot drift (Lane 2). Asserting
    full-string stability here would import two other lanes' open defects into this lane's
    guard, so the renumbering assertion is scoped to the nitrogen descriptors while
    rotation and reflection -- both of which ARE clean at whole-string level -- keep the
    strict form.
    """

    def _encode_variant(self, atoms, head):
        with tempfile.TemporaryDirectory() as d:
            return self._encode(_write_atoms(Path(d) / "probe.xyz", head, atoms))

    def test_descriptors_invariant_under_atom_renumbering(self):
        self._set(True)
        head, atoms = _read_atoms(RIFGUJ)
        base = sorted(_n_descriptors(self._encode_variant(atoms, head)))
        self.assertEqual(len(base), 3)
        rng = np.random.default_rng(11)
        for _ in range(4):
            perm = [int(i) for i in rng.permutation(len(atoms))]
            got = sorted(_n_descriptors(self._encode_variant([atoms[i] for i in perm], head)))
            self.assertEqual(got, base, "locked-N descriptors drifted under renumbering")

    def test_invariant_under_proper_rotation(self):
        self._set(True)
        head, atoms = _read_atoms(RIFGUJ)
        base = self._encode_variant(atoms, head)
        rng = np.random.default_rng(5)
        for _ in range(3):
            q, r = np.linalg.qr(rng.normal(size=(3, 3)))
            q = q @ np.diag(np.sign(np.diag(r)))
            if np.linalg.det(q) < 0:
                q[:, 0] *= -1
            self.assertEqual(self._encode_variant([(s, q @ c) for s, c in atoms], head), base)

    def test_flips_under_reflection(self):
        self._set(True)
        head, atoms = _read_atoms(RIFGUJ)
        base = self._encode_variant(atoms, head)
        mirrored = self._encode_variant(
            [(s, np.array([c[0], c[1], -c[2]])) for s, c in atoms], head
        )
        self.assertEqual(mirrored, _swap_n_tags(base))


class TestOverSensitivityGuard(_LeverBase):
    """Not every metal-bound nitrogen is a stereocentre, and most are not.

    An ammine (M,H,H,H) and a primary amine (M,H,H,R) have symmetry-equivalent hydrogens;
    emitting a configuration for them would be a fabricated descriptor the generator can
    never reproduce. Excluded by the four-distinct-symmetry-classes gate, so the check is
    made at both levels: no eligibility, and no change to the string.
    """

    ACHIRAL_N = (CISPLATIN, PT_AMMINE, PTCL2_EN)

    def _plan(self, path):
        mol, _ = get_tmc_mol(Path(path), 0, with_stereo=False)
        Chem.SanitizeMol(mol)
        CIPAssigner().assign_all(mol, diagnostics=False)
        return plan_locked_donors(mol)

    def test_ammine_and_primary_amine_are_not_eligible(self):
        for path in self.ACHIRAL_N:
            with self.subTest(fixture=Path(path).name):
                self.assertEqual(
                    self._plan(path).indices,
                    (),
                    "a metal-bound ammine / primary amine is NOT a stereocentre",
                )

    def test_ammine_and_primary_amine_strings_unchanged(self):
        for path in self.ACHIRAL_N + (BINAP,):
            with self.subTest(fixture=Path(path).name):
                self._set(False)
                off = self._encode(path)
                self._set(True)
                on = self._encode(path)
                self.assertEqual(off, on, "lever ON must not touch a non-stereogenic donor")
                self.assertNotIn("[N@", on)

    def test_pojjop_amine_is_eligible(self):
        """Positive control for the gate, so an over-tight gate cannot pass silently."""
        plan = self._plan(POJJOP)
        self.assertEqual(len(plan.indices), 1)


class TestPlanIsRenumberingInvariant(_LeverBase):
    """Eligibility itself must not depend on input numbering.

    The gate uses ``CanonicalRankAtoms(breakTies=False)``. ``breakTies=True`` was measured
    NOT invariant (2-11 distinct rank vectors over 20 renumberings), so this pins that the
    tie-free variant is used and stays used.
    """

    def test_eligible_count_stable_under_renumbering(self):
        head, atoms = _read_atoms(RIFGUJ)
        rng = np.random.default_rng(3)
        with tempfile.TemporaryDirectory() as d:
            counts = []
            for t in range(4):
                perm = [int(i) for i in rng.permutation(len(atoms))]
                p = _write_atoms(Path(d) / f"p{t}.xyz", head, [atoms[i] for i in perm])
                mol, _ = get_tmc_mol(Path(p), 0, with_stereo=False)
                Chem.SanitizeMol(mol)
                counts.append(len(plan_locked_donors(mol).indices))
        self.assertEqual(set(counts), {3}, f"eligibility drifted under renumbering: {counts}")


class TestRifgujRingCarbonsArePseudoAsymmetric(unittest.TestCase):
    """The guard that killed the naive P3-under-canonical-body fix. Keep it.

    RIFGUJ's three ring carbons are **pseudo-asymmetric**: they carry lowercase ``s``, a
    RELATIVE (all-cis) descriptor, and reflecting the molecule does NOT change them. Any change
    that makes the emitted string flip those carbons between a structure and its mirror is
    rewriting stereochemistry the geometry says is fixed -- which is exactly what happened when
    v0.4.6 tried to restore the locked-donor tag onto the reparsed fragment: setting a tag after
    the sanitize moved the canonical write order, and ``@``/``@@`` is a parity relative to that
    order.

    Runs in the SHIPPED configuration (all v0.4.5 defaults, locked-donor lever off), so it does
    not depend on either lever's state and cannot be pinned away.
    """

    def _cip(self, path):
        from rdkit.Chem import rdCIPLabeler

        mol, _ = get_tmc_mol(Path(path), 0, with_stereo=True)
        Chem.AssignStereochemistryFrom3D(mol)
        rdCIPLabeler.AssignCIPLabels(mol)
        return [
            (a.GetIdx(), a.GetProp("_CIPCode"))
            for a in mol.GetAtoms()
            if a.HasProp("_CIPCode") and a.GetSymbol() == "C"
        ]

    def test_geometry_says_the_ring_carbons_do_not_invert(self):
        with tempfile.TemporaryDirectory() as d:
            head, atoms = _read_atoms(RIFGUJ)
            mir = _write_atoms(
                Path(d) / "mirror.xyz",
                head,
                [(s, np.array([-c[0], c[1], c[2]])) for s, c in atoms],
            )
            base, mirror = self._cip(RIFGUJ), self._cip(mir)

        self.assertEqual(len(base), 3, f"expected 3 ring-carbon centres, got {base}")
        self.assertTrue(
            all(code.islower() for _, code in base),
            f"these must be pseudo-asymmetric (lowercase r/s), got {base}",
        )
        self.assertEqual(
            [c for _, c in base],
            [c for _, c in mirror],
            "reflection must NOT change a pseudo-asymmetric descriptor -- if this fails, either "
            "the fixture changed or something is now rewriting relative stereochemistry",
        )
