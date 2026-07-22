"""Guards for A3 (v0.4.3 vdw-acceptance): the generator must SEE van-der-Waals steric
clash, not just atomic fusion.

The release found 53% of generated structures carry >=1 vdW clash (vs 5% of real
crystals) because every acceptance gate judged only covalent-radius *fusion*. A3 adds a
vdW clash term that (1) shares one definition with the release quality metric
(``tools/structure_distortion_report.py``), (2) rejects a clashing conformer at embed
acceptance so the loop keeps searching for a clean one and the best-rejected fallback
surfaces the least-clashing candidate, and (3) ranks the pool + the geometry selection by
whole-complex clash.

Each test fails against pre-A3 code: the shared module did not exist, ``_finalize_positions``
accepted any non-fused conformer, and selection ranked on energy / coordination-sphere fit
alone.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np
from rdkit import Chem, RDLogger

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import clash, embed, om

RDLogger.DisableLog("rdApp.*")

CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"

_ADJ_RATIO = 1.4  # adj_ratio_criteria used by _finalize_positions
_ATOM_D = 0.5  # atom_d_criteria
_RATIO = 0.65  # ratio_criteria
_PT = Chem.GetPeriodicTable()


def _complex(oin):
    return om.get_om_from_modified_smiles(convert_parsed_to_msmiles(OINParser().parse(oin)))


def _finalize_args(mc):
    """Rebuild the (new_complex, radius_list, R, metal_index) get_embedding hands to
    _finalize_positions, so the criterion can be exercised on hand-crafted positions."""
    new_complex = mc.copy()
    atoms = new_complex.get_atom_list()
    radius_list = [a.get_radius() for a in atoms]
    n = len(radius_list)
    rr = np.repeat(np.array(radius_list), n).reshape((n, n))
    return new_complex, radius_list, rr + rr.T, new_complex.metal_index


def _finalize(mc, new_complex, radius_list, rmat, metal_index, positions):
    return embed._finalize_positions(
        positions.copy(),
        mc,
        new_complex,
        radius_list,
        rmat,
        metal_index,
        False,  # align: keep the hand-crafted coordinates as given
        _ATOM_D,
        _RATIO,
        _ADJ_RATIO,
    )


class TestClashHelperMatchesMetric(unittest.TestCase):
    """The gate's clash count must equal the release metric's on a real structure."""

    def test_helper_equals_geometry_metrics_on_broken_fixture(self):
        # gate == metric by construction, on a genuinely-clashing generated structure.
        sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
        from tools import structure_distortion_report as sdr

        path = os.path.join(os.path.dirname(__file__), "../fixtures/ticat3_generated_broken.xyz")
        lines = open(path).read().splitlines()
        nat = int(lines[0].split()[0])
        sym, xs = [], []
        for line in lines[2 : 2 + nat]:
            p = line.split()
            sym.append(p[0])
            xs.append([float(p[1]), float(p[2]), float(p[3])])
        xs = np.array(xs)
        z = [_PT.GetAtomicNumber(s) for s in sym]

        mine, _severe, _worst = clash.vdw_clash_count(xs, z)
        metric = sdr.geometry_metrics(sym, xs)["clash_vdw"]
        self.assertEqual(mine, metric, "gate clash count must match the release metric")
        self.assertGreater(metric, 0, "fixture is a known clashing structure")


class TestFinalizeVdwAcceptance(unittest.TestCase):
    """_finalize_positions accepts a clean conformer and rejects a topologically-valid
    but vdW-clashing one, scoring it into (-1, 0). The term is gated OFF by default (it
    loosens coordination on the current pool -- see clash.py), so these tests enable it."""

    def setUp(self):
        self._flag = clash.VDW_ACCEPTANCE_ENABLED
        clash.VDW_ACCEPTANCE_ENABLED = True
        self.mc = _complex(CISPLATIN_OIN)
        self.new_complex, self.radius_list, self.rmat, self.metal_index = _finalize_args(self.mc)
        clean = embed.get_embedding(self.mc, scale=1.0, option=0, align=True, seed=42)
        self.assertIsNotNone(clean, "cisplatin must embed to a clean geometry")
        self.clean = np.asarray(clean)
        self.atoms = self.new_complex.get_atom_list()
        self.elements = [a.get_element() for a in self.atoms]
        self.adj = np.asarray(self.mc.get_adj_matrix())

    def tearDown(self):
        clash.VDW_ACCEPTANCE_ENABLED = self._flag

    def test_disabled_by_default_accepts_clashing(self):
        # With the gate OFF (the shipped default), a vdW-clashing conformer that clears the
        # three covalent criteria is accepted exactly as pre-A3 -- proving the default path
        # is unchanged. (setUp turned it on; flip it off for this one assertion.)
        clash.VDW_ACCEPTANCE_ENABLED = False
        pos, _cl, _h = self._make_clash(2.0)
        self.assertGreater(
            clash.vdw_clash_count(pos, [_PT.GetAtomicNumber(e) for e in self.elements])[0], 0
        )
        good, candidate = _finalize(
            self.mc, self.new_complex, self.radius_list, self.rmat, self.metal_index, pos
        )
        self.assertIsNotNone(good, "gate OFF must accept the clashing conformer (pre-A3 behavior)")
        self.assertIsNone(candidate)

    def _make_clash(self, sep):
        """Move a Cl (bonded only to the metal -> diff-exempt) to ``sep`` A beyond an H,
        collinear with that H's N, so Cl-H clashes by vdW while every covalent bond and
        the perceived bond graph are unchanged."""
        cl = self.elements.index("Cl")
        h = self.elements.index("H")
        n = int(np.argmax(self.adj[h]))  # the H's bonded neighbour (its N)
        pos = self.clean.copy()
        u = pos[h] - pos[n]
        u = u / np.linalg.norm(u)
        pos[cl] = pos[h] + sep * u
        return pos, cl, h

    def test_clean_conformer_is_accepted(self):
        good, candidate = _finalize(
            self.mc, self.new_complex, self.radius_list, self.rmat, self.metal_index, self.clean
        )
        self.assertIsNotNone(good, "a clash-free conformer must still pass")
        self.assertIsNone(candidate)
        self.assertEqual(
            clash.vdw_clash_count(self.clean, [_PT.GetAtomicNumber(e) for e in self.elements])[0], 0
        )

    def test_clashing_conformer_is_rejected_and_scored(self):
        pos, cl, h = self._make_clash(2.0)
        z = [_PT.GetAtomicNumber(e) for e in self.elements]

        # Precondition: this geometry is exactly what pre-A3 code would ACCEPT --
        # it clears fusion, collapse and the adjacency-diff checks -- yet it clashes.
        dm = np.linalg.norm(pos[cl] - pos[h])
        self.assertGreater(dm, _ATOM_D, "not fused")
        rc = _PT.GetRcovalent(_PT.GetAtomicNumber("Cl")) + _PT.GetRcovalent(
            _PT.GetAtomicNumber("H")
        )
        self.assertGreater(dm / rc, _ADJ_RATIO, "not perceived as a new bond (diff stays 0)")
        clash_vdw, _s, _w = clash.vdw_clash_count(pos, z)
        self.assertGreater(clash_vdw, 0, "the crafted geometry really clashes by vdW")

        good, candidate = _finalize(
            self.mc, self.new_complex, self.radius_list, self.rmat, self.metal_index, pos
        )
        self.assertIsNone(good, "a vdW-clashing conformer must be rejected (accepted pre-A3)")
        self.assertIsNotNone(candidate)
        _returned_pos, score = candidate
        self.assertTrue(-1.0 < score < 0.0, f"clash score must sit in (-1, 0); got {score}")
        # Ranks above any topology-broken candidate (those score -diff <= -2), so the
        # unchanged best-rejected fallback prefers it.
        self.assertGreater(score, -2.0)

    def test_score_is_monotone_so_fallback_surfaces_least_clashing(self):
        z = [_PT.GetAtomicNumber(e) for e in self.elements]
        mild_pos, _, _ = self._make_clash(2.2)  # farther -> milder / fewer clashes
        severe_pos, _, _ = self._make_clash(1.6)  # closer -> worse overlap

        _g1, c_mild = _finalize(
            self.mc, self.new_complex, self.radius_list, self.rmat, self.metal_index, mild_pos
        )
        _g2, c_severe = _finalize(
            self.mc, self.new_complex, self.radius_list, self.rmat, self.metal_index, severe_pos
        )
        self.assertIsNotNone(c_mild)
        self.assertIsNotNone(c_severe)
        # Sanity: the closer placement is at least as clashy.
        self.assertLessEqual(
            clash.vdw_clash_count(mild_pos, z)[2], clash.vdw_clash_count(severe_pos, z)[2] + 1e-9
        )
        self.assertGreater(c_mild[1], c_severe[1], "milder clash must score higher (nearer 0)")
        # The fallback's argmax (embed.py) therefore returns the milder-clash candidate.
        self.assertEqual(max(c_mild[1], c_severe[1]), c_mild[1])


class TestPoolClashRanking(unittest.TestCase):
    """The shared clash key ranks least-clashing first -- the ordering used by the
    __init__ pool re-rank and _select_by_geometry."""

    def test_mol_clash_count_key_orders_least_clashing_first(self):
        class _Atom:
            def __init__(self, z, xyz):
                self._z, self._xyz = z, xyz

            def get_atomic_number(self):
                return self._z

            def get_coordinate(self):
                return np.array(self._xyz, dtype=float)

        class _Mol:
            def __init__(self, atom_list):
                self.atom_list = atom_list

        # Two carbons 2.3 A apart -> 1 vdW clash; 5 A apart -> clean.
        clashy = _Mol([_Atom(6, [0, 0, 0]), _Atom(6, [2.3, 0, 0])])
        clean = _Mol([_Atom(6, [0, 0, 0]), _Atom(6, [5.0, 0, 0])])
        self.assertEqual(clash.mol_clash_count(clean), 0)
        self.assertGreater(clash.mol_clash_count(clashy), 0)
        pool = [clashy, clean]
        pool.sort(key=clash.mol_clash_count)
        self.assertIs(pool[0], clean, "least-clashing conformer must rank first")


if __name__ == "__main__":
    unittest.main()
