"""Guard that the PuLP/CBC solver memo is a pure speed-up, not a behavior change.

MetalloGen's charge/bond-order perception (``compute_chg_and_bo`` ->
CBC subprocess) reads only the molecule's topology, yet is re-solved on the
identical graph across every conformer attempt. P2 memoizes it, keyed on
(atomic numbers + adjacency + charge + solver kwargs). Because the solve is a
pure function of that key, a cached answer is exactly what a fresh solve would
return -- so the generated XYZ must be *byte-identical* whether the cache is on
or off. These tests pin that contract:

* cache-on and cache-off produce byte-identical XYZ (the core correctness proof),
* a PuLP-bound golden actually exercises the cache (hits > 0 -- proof the
  redundant CBC subprocesses were collapsed), and
* generation stays deterministic across repeats with the cache on.

The fast, FF-only, small-pool config mirrors ``test_generation_determinism`` --
correctness of the memo is independent of pool width, so a 3-conformer pool
exercises the exact same solver path as the shipped default.
"""

import contextlib
import hashlib
import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import RDLogger

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generator3d.utils.compute_chg_and_bo_pulp import (
    clear_pulp_cache,
    pulp_cache_stats,
    set_pulp_cache_enabled,
)

RDLogger.DisableLog("rdApp.*")

# Small square-planar Pt(II) goldens reused across the suite; cisplatin is
# PuLP-bound (its wall-clock is dominated by the CBC solve), so it is the one
# guaranteed to register cache hits.
CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
TRANSPLATIN = "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"
CIS_PTCL2EN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}CCN{3}"
FLEXIBLE = "[Pt_SPL].C[C@H](O)CN{0}.[Cl]{1}.[Cl]{2}.N{3}"
GOLDEN_OINS = (CISPLATIN, TRANSPLATIN, CIS_PTCL2EN, FLEXIBLE)

# Small pool -> fast; the memo's correctness is independent of pool width.
_FAST = dict(engine="metallogen", optimizer="ff", ensemble_size=1, ff_params={"uff_pool_size": 3})


def _xyz_hash(oin, **kwargs):
    """SHA-256 of the generated XYZ for ``oin`` under the fast, FF-only config."""
    gen = OIN3DGenerator(**_FAST, **kwargs)
    with (
        open(os.devnull, "w") as dn,
        contextlib.redirect_stdout(dn),
        contextlib.redirect_stderr(dn),
    ):
        return hashlib.sha256(gen.generate(oin).xyz.encode()).hexdigest()


class TestPulpMemo(unittest.TestCase):
    def setUp(self):
        set_pulp_cache_enabled(True)
        clear_pulp_cache()

    def tearDown(self):
        # Never leave the global memo disabled for other tests in the process.
        set_pulp_cache_enabled(True)
        clear_pulp_cache()

    def test_cache_on_equals_cache_off(self):
        """The memo must not change the generated geometry for any golden."""
        for oin in GOLDEN_OINS:
            with self.subTest(oin=oin):
                set_pulp_cache_enabled(True)
                with_cache = _xyz_hash(oin)
                set_pulp_cache_enabled(False)
                without_cache = _xyz_hash(oin)
                set_pulp_cache_enabled(True)
                self.assertEqual(
                    with_cache,
                    without_cache,
                    "PuLP memo changed the output -- key is incomplete",
                )

    def test_cache_is_exercised(self):
        """A PuLP-bound golden must produce real cache hits (redundant solves collapsed)."""
        set_pulp_cache_enabled(True)
        clear_pulp_cache()
        _xyz_hash(CISPLATIN)
        stats = pulp_cache_stats()
        self.assertGreater(
            stats["hits"],
            0,
            f"cisplatin generation registered no PuLP cache hits: {stats}",
        )

    def test_deterministic_across_runs_with_cache(self):
        """Cache on, the same golden generates byte-identical XYZ across repeats."""
        set_pulp_cache_enabled(True)
        hashes = {_xyz_hash(CISPLATIN) for _ in range(3)}
        self.assertEqual(len(hashes), 1, "generation must stay deterministic with the memo on")


if __name__ == "__main__":
    unittest.main()
