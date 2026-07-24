"""Guards for SL1 (v0.4.4 acceptance-gate): generate-until-key-exact early-exit.

The default generator builds a pool of conformers then selects one by geometry classification;
it never checks the chosen conformer re-encodes to the same OIN. SL1 adds an opt-in
``early_exit`` (``OIN_EARLY_EXIT`` / ``ff_params["early_exit"]``): the attempt loop stops the
moment an embedded conformer INDEPENDENTLY re-encodes to the requested OIN's fac/mer key
(``accept_fn`` in ``generate_3d_structures``), and the adapter's accept-first pass returns it.

Each test fails against pre-SL1 code: ``generate_3d_structures`` had no ``accept_fn`` param and
``_select_by_geometry`` no ``early_exit`` param, so early-exit could not stop the pool build or
short-circuit selection.

These generate real 3D structures (integration-weight, but placed with the unit guards per the
SL1 handoff). Off-flag byte-identity to pristine is additionally covered by the determinism
goldens in tests/integration/test_roundtrip_smoke.py and tests/unit/test_regression_stability.py,
which run ``accept_fn=None`` unchanged.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import embed as embed_mod
from oinsmiles.generator3d import generate_3d_structures
from oinsmiles.oin.compare import normalize_oin_for_comparison, winding_canonical_key
from oinsmiles.utils.xyz2mol import get_oin_string

SEED = 42

CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
TRANSPLATIN = "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"
FERROCENE = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"


def _roundtrip_key(oin_string):
    return winding_canonical_key(normalize_oin_for_comparison(oin_string.strip()))


def _reencode(gen_mol, xyz):
    lines = xyz.splitlines()
    natoms = int(lines[0].strip())
    coords = np.array([[float(x) for x in lines[i].split()[1:4]] for i in range(2, 2 + natoms)])
    return get_oin_string(gen_mol, coords)


class _CountEmbeds:
    """Context manager that counts ``embed.get_embedding`` / ``get_embeddings_batch`` calls.

    ``generate_3d_structures`` references the ``embed`` module (``from . import embed``), so
    patching the module attribute counts calls made through it, serial or batched.
    """

    def __enter__(self):
        self.count = 0
        self._orig_single = embed_mod.get_embedding
        self._orig_batch = embed_mod.get_embeddings_batch
        outer = self

        def _single(*a, **k):
            outer.count += 1
            return outer._orig_single(*a, **k)

        def _batch(*a, **k):
            outer.count += int(k.get("num_confs", 1) or 1)
            return outer._orig_batch(*a, **k)

        embed_mod.get_embedding = _single
        embed_mod.get_embeddings_batch = _batch
        return self

    def __exit__(self, *exc):
        embed_mod.get_embedding = self._orig_single
        embed_mod.get_embeddings_batch = self._orig_batch
        return False


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


class TestPoolEarlyExit(unittest.TestCase):
    """The ``accept_fn`` hook stops the pool build early at the engine level."""

    def test_accept_fn_stops_building_and_returns_one(self):
        msmiles = _msmiles(CISPLATIN)

        with _CountEmbeds() as off:
            mols_off = generate_3d_structures(
                msmiles, num_conformers=5, uff_pool_size=10, seed=SEED
            )
        with _CountEmbeds() as on:
            mols_on = generate_3d_structures(
                msmiles,
                num_conformers=5,
                uff_pool_size=10,
                seed=SEED,
                accept_fn=lambda _m: True,  # accept the very first cleaned conformer
            )

        self.assertEqual(len(mols_on), 1, "early-exit returns exactly the matched conformer")
        self.assertGreater(len(mols_off), 1, "off-flag builds a multi-conformer pool")
        self.assertLess(on.count, off.count, "early-exit embeds fewer conformers than full-pool")

    def test_accept_fn_none_is_unchanged(self):
        # accept_fn=None (default) must build the pool exactly as before.
        msmiles = _msmiles(CISPLATIN)
        with _CountEmbeds() as a:
            mols_a = generate_3d_structures(msmiles, num_conformers=5, uff_pool_size=10, seed=SEED)
        with _CountEmbeds() as b:
            mols_b = generate_3d_structures(
                msmiles, num_conformers=5, uff_pool_size=10, seed=SEED, accept_fn=None
            )
        self.assertEqual(a.count, b.count)
        self.assertEqual(len(mols_a), len(mols_b))


class TestEarlyExitKeyExact(unittest.TestCase):
    """OIN_EARLY_EXIT returns a key-exact conformer, in fewer embed attempts."""

    GOLDENS = {"cisplatin": CISPLATIN, "transplatin": TRANSPLATIN, "ferrocene": FERROCENE}

    def test_goldens_key_exact_and_fewer_attempts(self):
        total_off = 0
        total_on = 0
        for name, oin in self.GOLDENS.items():
            with self.subTest(fixture=name):
                with _CountEmbeds() as off:
                    # off-flag run: we only need its embed count for the comparison.
                    # early_exit is default-ON since v0.4.4, so opt OUT explicitly here to
                    # measure the full-pool baseline.
                    OIN3DGenerator(
                        optimizer="ff", seed=SEED, ff_params={"early_exit": False}
                    ).generate(oin)
                with _CountEmbeds() as on:
                    gen_on = OIN3DGenerator(
                        optimizer="ff", seed=SEED, ff_params={"early_exit": True}
                    ).generate(oin)

                self.assertIsNotNone(gen_on.mol, f"{name}: early-exit returned no bonded mol")
                key_in = _roundtrip_key(oin)
                key_on = _roundtrip_key(_reencode(gen_on.mol, gen_on.xyz))
                self.assertEqual(key_in, key_on, f"{name}: early-exit conformer is not key-exact")
                # Never MORE attempts than the full-pool path.
                self.assertLessEqual(
                    on.count, off.count, f"{name}: early-exit used more embeds than full-pool"
                )
                total_off += off.count
                total_on += on.count

        # Aggregate: early-exit genuinely saves embeds across the goldens.
        self.assertLess(total_on, total_off, "early-exit saves embeds in aggregate over goldens")

    def test_early_exit_on_by_default(self):
        # v0.4.4 promote A/B flipped early-exit to default-ON (+15.8pt byte-exact, 0
        # regressions, ~5x faster). Guard: a DEFAULT run (no ff_params) must early-exit --
        # i.e. use fewer embeds in aggregate than an explicit early_exit=False full-pool run.
        # Fails if the default is silently reverted to OFF.
        total_default = 0
        total_off = 0
        for oin in self.GOLDENS.values():
            with _CountEmbeds() as dflt:
                OIN3DGenerator(optimizer="ff", seed=SEED).generate(oin)
            with _CountEmbeds() as off:
                OIN3DGenerator(optimizer="ff", seed=SEED, ff_params={"early_exit": False}).generate(
                    oin
                )
            total_default += dflt.count
            total_off += off.count
        self.assertLess(
            total_default, total_off, "early_exit must be ON by default (v0.4.4 promotion)"
        )


if __name__ == "__main__":
    unittest.main()
