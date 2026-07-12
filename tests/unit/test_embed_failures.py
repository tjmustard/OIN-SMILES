"""Guards for P9 (v0.4.0-perf embed-failures): the primary ETKDG embed no longer
returns -1 on the loose-scale OCT tris-chelate combos.

Root cause (measured with params.trackFailures / GetFailureCounts): with
useBasicKnowledge=False the distance-geometry embed carries no knowledge terms, so
`useRandomCoords=True` hands the first minimization an infeasible random start and it
fails FIRST_MINIMIZATION on all 100 iterations, returning -1 after paying the full
maxIterations run. On fac-Ir(ppy)3 the 4 seed-42 failing combos are all LOOSE scales
(1.1/1.2 x option 0/1) -- a looser cmap over-stretches the fused-aromatic chelate bite.

P9 sets `useRandomCoords=False` (seed the embed from the metric-matrix eigen-
decomposition, RDKit's default) so it starts near-feasible: the failures flip to
successes AND the embed gets ~3x cheaper. These tests pin both the flip and that the
seeded embed stays reproducible under the new params.

Each test fails against pre-P9 code (useRandomCoords=True), where the loose-scale
combos return -1 / None.
"""

import itertools
import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np
from rdkit import RDLogger

from oinsmiles.generation.metallogen_adapter import convert_parsed_to_msmiles
from oinsmiles.generation.oin_parser import OINParser
from oinsmiles.generator3d import embed, om

RDLogger.DisableLog("rdApp.*")

# fac-Ir(ppy)3: CN6 octahedral tris-chelate, the exemplar hard-embed OCT complex.
IR_OIN = "[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"
SCALES = [0.8, 0.9, 1.0, 1.1, 1.2]
OPTIONS = [0, 1, 2]
# The (scale, option) combos that return -1 on pre-P9 code (useRandomCoords=True),
# measured at seed 42 on fac-Ir(ppy)3.
LOOSE_FAILING_COMBOS = [(1.1, 0), (1.1, 1), (1.2, 0), (1.2, 1)]


def _complex(oin):
    return om.get_om_from_modified_smiles(convert_parsed_to_msmiles(OINParser().parse(oin)))


class TestPrimaryEmbedFailures(unittest.TestCase):
    def test_loose_scale_failures_are_rescued_by_the_retry(self):
        """P9: every (scale, option) combo yields a conformer on fac-Ir(ppy)3 at seed 42.

        The primary random-coords embed still returns -1 on the loose-scale combos (that
        is unchanged and expected); what P9 adds is a metric-matrix retry
        (useRandomCoords=False) that rescues each one, so no combo ends up without a
        conformer. This asserts (a) zero outright failures, and (b) that the rescue
        actually came from the retry -- at least one primary -1 occurred and every such
        -1 was followed by a successful useRandomCoords=False retry. The primary vs retry
        calls are told apart by params.useRandomCoords at call time.
        """
        mc = _complex(IR_OIN)
        primary_fail = {"n": 0}
        retry_success = {"n": 0}
        orig = embed.AllChem.EmbedMolecule

        def spy(*a, **k):
            params = a[1]
            is_retry = not params.useRandomCoords
            rc = orig(*a, **k)
            if not is_retry and rc == -1:
                primary_fail["n"] += 1
            if is_retry and rc != -1:
                retry_success["n"] += 1
            return rc

        embed.AllChem.EmbedMolecule = spy
        none_count = 0
        try:
            for scale, option in itertools.product(SCALES, OPTIONS):
                pos = embed.get_embedding(mc.copy(), scale, option, align=True, seed=42)
                none_count += int(pos is None)
        finally:
            embed.AllChem.EmbedMolecule = orig

        self.assertEqual(none_count, 0, "every combo must yield a conformer after the P9 retry")
        self.assertGreater(
            primary_fail["n"], 0, "expected the loose-scale primary embeds to still return -1"
        )
        self.assertGreaterEqual(
            retry_success["n"],
            primary_fail["n"],
            "every primary -1 must be rescued by a successful useRandomCoords=False retry",
        )

    def test_loose_scale_combos_now_return_a_conformer(self):
        """The specific pre-P9 failing combos now yield a valid embedding (not None)."""
        mc = _complex(IR_OIN)
        for scale, option in LOOSE_FAILING_COMBOS:
            with self.subTest(scale=scale, option=option):
                pos = embed.get_embedding(mc.copy(), scale, option, align=True, seed=42)
                self.assertIsNotNone(
                    pos, f"loose combo (scale={scale}, option={option}) still fails to embed"
                )
                self.assertEqual(len(pos), mc.num_atom, "positions must cover every atom")

    def test_embed_is_deterministic_under_new_params(self):
        """useRandomCoords=False must stay reproducible: same seed -> identical coords."""
        mc = _complex(IR_OIN)
        a = embed.get_embedding(mc.copy(), 1.1, 0, align=True, seed=42)
        b = embed.get_embedding(mc.copy(), 1.1, 0, align=True, seed=42)
        self.assertIsNotNone(a)
        np.testing.assert_array_equal(
            np.asarray(a), np.asarray(b), "same seed must give byte-identical embed coords"
        )


if __name__ == "__main__":
    unittest.main()
