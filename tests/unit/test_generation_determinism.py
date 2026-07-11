"""Guard that MetalloGen generation is reproducible via an explicit ``seed``.

The dummy-metal ETKDG embed used to draw its ``params.randomSeed`` from an
unseeded ``random.randint``, so the same OIN produced a different 3D structure on
every run. S6 made ``seed=42`` the deterministic default at the embed layer; this
phase (P1) threads a caller-chosen ``seed`` the rest of the way up
(``OIN3DGenerator`` -> ``MetalloGenAdapter`` -> ``generate_3d_structures``) and
exposes it as ``oin-smiles oin2xyz --seed``. These tests pin that contract:

* the same seed yields byte-identical XYZ across runs,
* the default (no explicit seed) is deterministic and equals ``seed=42``,
* two different seeds yield different structures (the seed reaches ETKDG), and
* the per-attempt seed offset keeps the conformer pool diverse (it does not
  collapse to N copies of one conformer).

Most cases use a deliberately small UFF pool (``uff_pool_size``/``ensemble_size``)
purely to keep the guard fast: determinism is a property of the seed, not of the
pool width, so a 3-conformer pool exercises the exact same seed->embed thread as
the shipped default. The default-pool, 5-run byte-identity on all four goldens
was verified manually while landing P1.
"""

import contextlib
import hashlib
import os
import subprocess
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import RDLogger

from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.metallogen_adapter import convert_oin_to_msmiles
from oinsmiles.generator3d import calculate_heavy_atom_rmsd, generate_3d_structures

RDLogger.DisableLog("rdApp.*")

# Four small, fast square-planar Pt(II) goldens reused across the suite. The last
# has a flexible amino-alcohol arm, so different seeds give visibly different
# geometry (rigid cisplatin can converge to the same structure regardless).
CISPLATIN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
TRANSPLATIN = "[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}"
CIS_PTCL2EN = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}CCN{3}"
FLEXIBLE = "[Pt_SPL].C[C@H](O)CN{0}.[Cl]{1}.[Cl]{2}.N{3}"
GOLDEN_OINS = (CISPLATIN, TRANSPLATIN, CIS_PTCL2EN, FLEXIBLE)

# Small pool -> fast; determinism is independent of pool width (see module docstring).
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


class TestSeedDeterminism(unittest.TestCase):
    def test_same_seed_byte_identical(self):
        """Each golden OIN generates byte-identical XYZ across runs at a fixed seed."""
        for oin in GOLDEN_OINS:
            with self.subTest(oin=oin):
                hashes = {_xyz_hash(oin, seed=42) for _ in range(3)}
                self.assertEqual(len(hashes), 1, "seed=42 must be byte-reproducible")

    def test_default_seed_is_deterministic_and_is_42(self):
        """No explicit seed is deterministic (default 42), not random per run."""
        default_hashes = {_xyz_hash(CISPLATIN) for _ in range(2)}
        self.assertEqual(len(default_hashes), 1, "the default seed must be deterministic")
        self.assertEqual(
            default_hashes.pop(),
            _xyz_hash(CISPLATIN, seed=42),
            "the default must equal an explicit seed=42",
        )

    def test_different_seeds_differ(self):
        """Two seeds give two structures -- proof the seed reaches the ETKDG embed."""
        self.assertNotEqual(
            _xyz_hash(FLEXIBLE, seed=7),
            _xyz_hash(FLEXIBLE, seed=42),
            "different seeds must produce different geometry",
        )

    def test_pool_diversity_under_seed(self):
        """The per-attempt seed offset keeps the pool diverse (no collapse to copies)."""
        msmiles = convert_oin_to_msmiles(FLEXIBLE)
        with (
            open(os.devnull, "w") as dn,
            contextlib.redirect_stdout(dn),
            contextlib.redirect_stderr(dn),
        ):
            pool = generate_3d_structures(
                msmiles, num_conformers=8, optimizer=None, uff_pool_size=8, seed=42
            )
        self.assertGreaterEqual(len(pool), 2, "seeded pool collapsed to a single conformer")
        max_rmsd = max(calculate_heavy_atom_rmsd(pool[0], m) for m in pool[1:])
        self.assertGreater(max_rmsd, 0.5, "seeded pool members are not distinct conformers")


class TestCliSeed(unittest.TestCase):
    def test_cli_threads_seed_reproducibly(self):
        """oin2xyz --seed threads through to a reproducible XYZ block."""

        def run_cli():
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "oinsmiles.cli",
                    "oin2xyz",
                    CISPLATIN,
                    "--optimizer",
                    "ff",
                    "--seed",
                    "7",
                ],
                capture_output=True,
                text=True,
                timeout=180,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr[-2000:])
            return proc.stdout

        self.assertEqual(run_cli(), run_cli(), "oin2xyz --seed must be reproducible")


if __name__ == "__main__":
    unittest.main()
