"""Guards for P3 (v0.4.0-perf etkdg-embed): the ETKDG return-code check, the removal
of the dead rebuild-retry ladder, and the opt-in C++-parallel batched embed.

Two independent wins, each pinned here:

* ``EmbedMolecule`` returns ``-1`` on failure; it does not raise. The serial embed
  now checks that return code directly and, on a plain ``-1``, does NOT run the
  rebuild retry -- that retry re-embeds the identical PuLP mol with the identical
  seeded ``params`` and reproduces the same ``-1`` (measured 8/8 dead on
  fac-Ir(ppy)3). Its ``continue`` also used to suppress the ``if not haptic_exist:
  break``, forcing 3 further redundant primary embeds per failing non-haptic combo.
  Removing it is byte-identical (proven by the golden A/B) and ~3.5x faster on Ir.
  (P9 later re-introduced a *different*, materially-distinct retry on -1 -- from
  metric-matrix initial coords, ``useRandomCoords=False`` -- which flips the OCT
  loose-scale failures instead of reproducing the same -1; that is why the all-(-1)
  count below is 2x per alt mol, not 1x. It is deterministic and fires only on -1, so
  every succeeding embed keeps its exact geometry. See test_embed_failures.py.)

* ``get_embeddings_batch`` + ``embed_num_threads`` add an opt-in batched embed via
  ``EmbedMultipleConfs`` (``num_threads=0`` = all cores). It is gated behind the
  knob so the default (``num_threads=1``) stays the serial, byte-identical path.

Each test fails against the pre-P3 code (the return code was never checked; the
rebuild ran on every ``-1``; there was no batched path or ``num_threads`` knob).
"""

import hashlib
import os
import sys
import unittest
from unittest import mock

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import RDLogger

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import embed, generate_3d_structures, om

RDLogger.DisableLog("rdApp.*")

CISPLATIN_OIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"  # CN4 square planar, non-haptic
FERROCENE_OIN = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"  # haptic
IR_OIN = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"


def _msmiles(oin):
    return convert_parsed_to_msmiles(OINParser().parse(oin))


def _complex(oin):
    return om.get_om_from_modified_smiles(_msmiles(oin))


def _xyz_sha(oin, **ff):
    mols = generate_3d_structures(_msmiles(oin), num_conformers=1, optimizer=None, **ff)
    from oinsmiles.generator3d import get_xyz_string

    return hashlib.sha256(get_xyz_string(mols[0]).encode()).hexdigest()


class TestReturnCodeAndDeadRetry(unittest.TestCase):
    """The ``-1`` return code drives failure, and the rebuild retry is gone for it."""

    def test_minus_one_return_skips_rebuild_and_haptic_reloop(self):
        # With every EmbedMolecule returning -1 (never raising), get_embedding must
        # give up quickly. Pre-P3 this was ~8x per alt mol (4 haptic scales x {primary,
        # rebuild}). P3 cut it to exactly ONE primary embed per alt mol. P9 then added
        # ONE deterministic, materially-different retry from metric-matrix initial
        # coords (useRandomCoords=False) that fires only on a plain -1 -- so the count
        # is now exactly TWO per alt mol (primary + that retry). What stays removed:
        # the metal-swap rebuild ladder (raised-embed path only) and the 4x
        # scales_for_haptic re-loop on this non-haptic complex. Since the mock returns
        # -1 regardless of params, the P9 retry also fails here, so no valid position.
        mc = _complex(CISPLATIN_OIN)
        n_alt = len(embed.get_alternative_molecule(mc.copy(), 0)[0])
        calls = {"n": 0}
        orig = embed.AllChem.EmbedMolecule

        def always_fail(*a, **k):
            calls["n"] += 1
            return -1  # -1, NOT an exception

        embed.AllChem.EmbedMolecule = always_fail
        try:
            positions = embed.get_embedding(mc, scale=1.0, option=0, align=True, seed=42)
        finally:
            embed.AllChem.EmbedMolecule = orig

        self.assertIsNone(positions, "every embed returned -1, so there is no valid position")
        self.assertEqual(
            calls["n"],
            2 * n_alt,
            "a -1 must cost one primary embed + one P9 useRandomCoords=False retry per "
            "alt mol: no metal-swap rebuild, no haptic re-loop",
        )

    def test_minus_one_not_treated_as_success(self):
        # A -1 return with NO exception must still be detected as a failure. If the
        # code relied only on GetConformer() raising, a stubbed embed that returns -1
        # while leaving a stale/again conformer could leak through; the rc check
        # closes that. Here the whole generation must yield nothing.
        with mock.patch.object(embed.AllChem, "EmbedMolecule", lambda *a, **k: -1):
            res = generate_3d_structures(
                _msmiles(CISPLATIN_OIN), ff_params={"max_attempts": 5}, optimizer=None
            )
        self.assertEqual(res, [], "no embed produced a conformer, so the pool must be empty")


class TestNumThreadsDefaultIsSerial(unittest.TestCase):
    """Default num_threads=1 keeps the serial path; it is byte-identical to explicit 1."""

    def test_default_does_not_take_the_batched_path(self):
        seen = {"batched": 0}
        real_batch = embed.get_embeddings_batch

        def spy(*a, **k):
            seen["batched"] += 1
            return real_batch(*a, **k)

        with mock.patch.object(embed, "get_embeddings_batch", spy):
            generate_3d_structures(_msmiles(CISPLATIN_OIN), optimizer=None)
        self.assertEqual(seen["batched"], 0, "the default must not call the batched embed")

    def test_num_threads_zero_takes_the_batched_path(self):
        seen = {"batched": 0}
        real_batch = embed.get_embeddings_batch

        def spy(*a, **k):
            seen["batched"] += 1
            return real_batch(*a, **k)

        with mock.patch.object(embed, "get_embeddings_batch", spy):
            generate_3d_structures(
                _msmiles(CISPLATIN_OIN), optimizer=None, ff_params={"embed_num_threads": 0}
            )
        self.assertGreater(seen["batched"], 0, "embed_num_threads=0 must use the batched embed")

    def test_default_equals_explicit_one_byte_for_byte(self):
        self.assertEqual(
            _xyz_sha(CISPLATIN_OIN),
            _xyz_sha(CISPLATIN_OIN, ff_params={"embed_num_threads": 1}),
            "num_threads=1 must be the default and produce identical geometry",
        )

    def test_embed_num_threads_is_stripped_from_cleaner_params(self):
        # embed_num_threads is a generation knob, not a TMCOptimizer arg; it must be
        # consumed by the fill loop and not forwarded into TMCOptimizer(**...).
        res = generate_3d_structures(
            _msmiles(CISPLATIN_OIN), optimizer=None, ff_params={"embed_num_threads": 0}
        )
        self.assertTrue(res, "generation with embed_num_threads set must still succeed")


class TestBatchedEmbed(unittest.TestCase):
    """get_embeddings_batch: haptic fallback signal, valid output, determinism, diversity."""

    def test_returns_none_for_haptic(self):
        batch = embed.get_embeddings_batch(
            _complex(FERROCENE_OIN), scale=1.0, option=0, num_confs=4, num_threads=0, seed=42
        )
        self.assertIsNone(batch, "haptic complexes must signal serial fallback (None)")

    def test_returns_positions_for_non_haptic(self):
        batch = embed.get_embeddings_batch(
            _complex(CISPLATIN_OIN), scale=1.0, option=0, num_confs=4, num_threads=0, seed=42
        )
        self.assertIsInstance(batch, list)
        self.assertGreaterEqual(len(batch), 1, "a feasible non-haptic combo must yield positions")
        mc = _complex(CISPLATIN_OIN)
        for positions in batch:
            self.assertEqual(len(positions), mc.num_atom, "positions must be trimmed to num_atom")

    def test_num_threads_zero_is_deterministic(self):
        self.assertEqual(
            _xyz_sha(CISPLATIN_OIN, ff_params={"embed_num_threads": 0}),
            _xyz_sha(CISPLATIN_OIN, ff_params={"embed_num_threads": 0}),
            "the same seed must give byte-identical XYZ even with parallel embedding",
        )

    def test_batched_pool_stays_diverse(self):
        # A faster generator that collapses the pool to one conformer would break
        # _select_by_geometry. Ir(ppy)3 must still yield >= 2 distinct conformers.
        from oinsmiles.generator3d import calculate_heavy_atom_rmsd

        mols = generate_3d_structures(
            _msmiles(IR_OIN),
            num_conformers=6,
            optimizer=None,
            uff_pool_size=6,
            ff_params={"embed_num_threads": 0},
        )
        kept = []
        for m in mols:
            if all(calculate_heavy_atom_rmsd(m, k) >= 0.5 for k in kept):
                kept.append(m)
        self.assertGreaterEqual(len(kept), 2, "batched Ir(ppy)3 pool must keep >= 2 distinct confs")


if __name__ == "__main__":
    unittest.main()
